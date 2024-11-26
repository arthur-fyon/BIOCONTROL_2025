#=
This file contains large simulations scripts
=#

using DifferentialEquations, ProgressMeter

include("STG_kinetics.jl") # Include STG kinetics of gating variables
include("STG_models.jl") # Include STG model
include("STG_utils.jl") # Include some utils functions
include("STG_gs_derivatives.jl") # Include X_inf derivatives
include("STG_DIC.jl") # Include the DIC and compensation algorithm
include("STG_neuromodulation.jl") # Include the neuromodulation cells functions

# Moving average function
moving_average(vs, n, padding) = [sum(vs[i:(i+n-1)])/n for i in 1:padding:(length(vs)-(n-1))];

## STG model from Liu 1998 - current-clamp mode
function simulate_STG(g_all, Iapp, tau_Na, tau_g, Ca_tgt, C)
    ## Simulation of the model in current-clamp mode
    # Extracting the maximal ion channel conductances
    (gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak) = g_all

    # Computing homeostasis time constants
    tau_CaT = tau_Na * gNa / gCaT
    tau_KCa = tau_Na * gNa / gKCa
    tau_Kd = tau_Na * gNa / gKd
    tau_H = tau_Na * gNa / gH
    tau_leak = tau_Na * gNa / gleak
    tau_CaS = tau_Na * gNa / gCaS
    tau_A = tau_Na * gNa / gA

    # Parameter vector for simulations
    p = [Iapp, tau_Na, tau_CaT, tau_CaS, tau_A, tau_KCa, tau_Kd, tau_H, tau_leak, tau_g, Ca_tgt, C]

     # Initial conditions
     V0  = -70.
     Ca0 = 0.5
     x0  = [V0, mNa_inf(V0), hNa_inf(V0), mCaT_inf(V0), hCaT_inf(V0), mCaS_inf(V0), hCaS_inf(V0),
          mA_inf(V0), hA_inf(V0), mKCa_inf(V0, Ca0), mKd_inf(V0), mH_inf(V0), Ca0,
          gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak,
          gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak]

    # Simulation
    prob = ODEProblem(STG_homeo_leak_ODE, x0, tspan_half, p) # Describing the problem
    sol = solve(prob) # Solving the problem

    return sol
end

## STG model from Liu 1998 - current-clamp mode
function simulate_STG_crash(g_all_crash, Iapp, tau_Na, tau_g, Ca_tgt, C, Vth)
    ## Simulation of the model in current-clamp mode
    # Extracting the maximal ion channel conductances
    (gNa_crash, gCaT_crash, gCaS_crash, gA_crash, gKCa_crash, gKd_crash, gH_crash, gleak_crash) = g_all_crash

    (gCaS_crash_nmod, gA_crash_nmod) = DICs_gmax_neuromodCaSA(gNa_crash, gCaT_crash,
        gKd_crash, gKCa_crash, gH_crash, gleak_crash, -8., 4., Vth)

    # Computing homeostasis time constants
    tau_CaT_crash = tau_Na * gNa_crash / gCaT_crash
    tau_KCa_crash = tau_Na * gNa_crash / gKCa_crash
    tau_Kd_crash = tau_Na * gNa_crash / gKd_crash
    tau_H_crash = tau_Na * gNa_crash / gH_crash
    tau_leak_crash = tau_Na * gNa_crash / gleak_crash
    tau_CaS_crash = tau_Na * gNa_crash / gCaS_crash
    tau_A_crash = tau_Na * gNa_crash / gA_crash

    # Parameter vector for simulations
    p_crash = [Iapp, tau_Na, tau_CaT_crash, tau_CaS_crash, tau_A_crash, tau_KCa_crash,
               tau_Kd_crash, tau_H_crash, tau_leak_crash, tau_g, Ca_tgt, C]

     # Initial conditions
     V0  = -70.
     Ca0 = 0.5
     x0  = [V0, mNa_inf(V0), hNa_inf(V0), mCaT_inf(V0), hCaT_inf(V0), mCaS_inf(V0), hCaS_inf(V0),
          mA_inf(V0), hA_inf(V0), mKCa_inf(V0, Ca0), mKd_inf(V0), mH_inf(V0), Ca0,
          gNa_crash, gCaT_crash, gCaS_crash_nmod, gA_crash_nmod, gKCa_crash, gKd_crash, gH_crash, gleak_crash,
          gNa_crash, gCaT_crash, gCaS_crash_nmod, gA_crash_nmod, gKCa_crash, gKd_crash, gH_crash, gleak_crash]

    # Simulation
    prob = ODEProblem(STG_homeo_leak_ODE, x0, tspan_half, p_crash) # Describing the problem
    sol = solve(prob) # Solving the problem

    return sol
end

## STG model from Liu 1998 - current-clamp mode
function simulate_STG_population(g_all_init, Iapp, tau_Na, tau_g, Ca_tgt,
                                 C, ICs_th_init, tt_short)
    # Initializing some variables
    burstiness_after = zeros(ncells)
    gNa_matrix = zeros(ncells, length(tt[2:5000:end]))
    gCaT_matrix = zero(gNa_matrix)
    gCaS_matrix = zero(gNa_matrix)
    gA_matrix = zero(gNa_matrix)
    gKCa_matrix = zero(gNa_matrix)
    gKd_matrix = zero(gNa_matrix)
    gH_matrix = zero(gNa_matrix)
    gleak_matrix = zero(gNa_matrix)

    window_size_s = 5
    padding_s = 0.2
    window_size_ = Int64(round(window_size_s*1000*length(tt)/(Tfinal)))
    padding_ = Int64(round(padding_s*1000*length(tt)/(Tfinal)))
    tt_moving_average_ = window_size_s/2 : padding_s : Tfinal/1000 - window_size_s/2

    Ca_ma_matrix = zeros(ncells, length(tt_moving_average_))

    @showprogress "Computing..." for i = 1 : ncells
        # Simulate to initial blockade conditions
        sol = simulate_STG(g_all_init[i, :], Iapp, tau_Na, tau_g, Ca_tgt, C)

        # Retrieving variables
        x = sol(tt_short[2:5000:end])

        # Simulate to crash conditions
        g_all_crash = [x[14, end], x[15, end], x[16, end], x[17, end], x[18, end], x[19, end], x[20, end], x[21, end]]
        sol_crash = simulate_STG_crash(g_all_crash, Iapp, tau_Na, tau_g, Ca_tgt, C, ICs_th_init[:, 1][i])

        # Computing both burstiness
        V_burstiness = sol_crash(tt_short)[1, :]
        i_after = findall(tt_short .> 10000 .&& tt_short .< 15000)
        burstiness_after[i], _, _, _ = extract_burstiness(V_burstiness[i_after], tt_short[i_after])

        # Retrieving variables
        x_crash = sol_crash(tt_short[2:5000:end])

        # Pushing variables into matrices
        gNa_matrix[i, :] = [x[14, :]; x_crash[14, :]]
        gCaT_matrix[i, :] = [x[15, :]; x_crash[15, :]]
        gCaS_matrix[i, :] = [x[16, :]; x_crash[16, :]]
        gA_matrix[i, :] = [x[17, :]; x_crash[17, :]]
        gKCa_matrix[i, :] = [x[18, :]; x_crash[18, :]]
        gKd_matrix[i, :] = [x[19, :]; x_crash[19, :]]
        gH_matrix[i, :] = [x[20, :]; x_crash[20, :]]
        gleak_matrix[i, :] = [x[21, :]; x_crash[21, :]]

        Ca_all = [sol(tt_short)[13, :]; sol_crash(tt_short)[13, :]]
        Ca_ma_matrix[i, :] = moving_average(Ca_all, window_size_, padding_)
    end

    return gNa_matrix, gCaT_matrix, gCaS_matrix, gA_matrix, gKCa_matrix,
           gKd_matrix, gH_matrix, gleak_matrix, Ca_ma_matrix, burstiness_after
end
