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
function simulate_STG(g_all, Iapp, tau_Na, tau_g, Ca_tgt, C, α, β, Kp, Ki, Kt,
                      gsth_sim, guth_sim, Vth, u_maxCaS, u_maxA)
    ## Simulation of the model in current-clamp mode
    # Extracting the maximal ion channel conductances
    (gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak) = g_all

    # Computing homeostasis time constants
    tau_CaT = tau_Na * gNa / gCaT
    tau_KCa = tau_Na * gNa / gKCa
    tau_Kd = tau_Na * gNa / gKd
    tau_H = tau_Na * gNa / gH
    tau_leak = tau_Na * gNa / gleak

    # Parameter vector for simulations
    p = [Iapp, tau_Na, tau_CaT, tau_KCa, tau_Kd, tau_H, tau_leak, tau_g, Ca_tgt, C, α, β,
         Kp, Ki, Kt, gsth_sim, guth_sim, Vth, u_maxCaS, u_maxA]

     # Initial conditions
     V0  = -70.
     Ca0 = 0.5
     x0  = [V0, mNa_inf(V0), hNa_inf(V0), mCaT_inf(V0), hCaT_inf(V0), mCaS_inf(V0),
            hCaS_inf(V0), mA_inf(V0), hA_inf(V0), mKCa_inf(V0, Ca0), mKd_inf(V0), mH_inf(V0),
            Ca0, gNa, gCaT, gKCa, gKd, gH, gleak, gNa, gCaT, gKCa, gKd, gH, gleak,
            gCaS, gCaS, (β * gCaS) / Ki, gA, gA, (β * gA) / Ki]

    # Simulation
    prob = ODEProblem(simple_PI_homeo_leak_STG_ODE, x0, tspan, p) # Describing the problem
    sol = solve(prob) # Solving the problem

    return sol
end

## STG model from Liu 1998 - current-clamp mode
function simulate_STG_population(g_all_init, Iapp, tau_Na, tau_g, Ca_tgt,
                                 C, α, β, Kp, Ki, Kt, gsth_sim, guth_sim,
                                 ICs_th_init, u_maxCaS, u_maxA, tt)
    # Initializing some variables
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
        sol = simulate_STG(g_all_init[i, :], Iapp, tau_Na, tau_g, Ca_tgt, C, α, β, Kp, Ki, Kt,
                           gsth_sim, guth_sim, ICs_th_init[:, 1][i], u_maxCaS, u_maxA)

        # Retrieving variables
        x = sol(tt[2:5000:end])

        # Pushing variables into matrices
        gNa_matrix[i, :] = x[14, :]
        gCaT_matrix[i, :] = x[15, :]
        gCaS_matrix[i, :] = x[27, :]
        gA_matrix[i, :] = x[30, :]
        gKCa_matrix[i, :] = x[16, :]
        gKd_matrix[i, :] = x[17, :]
        gH_matrix[i, :] = x[18, :]
        gleak_matrix[i, :] = x[19, :]

        Ca_ma_matrix[i, :] = moving_average(sol(tt)[13, :], window_size_, padding_)
    end

    return gNa_matrix, gCaT_matrix, gCaS_matrix, gA_matrix, gKCa_matrix,
           gKd_matrix, gH_matrix, gleak_matrix, Ca_ma_matrix
end
