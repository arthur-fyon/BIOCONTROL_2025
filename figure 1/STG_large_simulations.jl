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
function simulate_STG_first_half(g_all, Iapp, tau_Na, tau_CaT, tau_CaS, tau_A, tau_KCa, tau_Kd, tau_H, tau_leak, Ca_tgt, C)
    ## Simulation of the model in current-clamp mode
    # Extracting the maximal ion channel conductances
    (gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak) = g_all

    # Parameter vector for simulations
    p = (Iapp, tau_Na, tau_CaT, tau_CaS, tau_A, tau_KCa, tau_Kd, tau_H, tau_leak, tau_g, Ca_tgt, C)

     # Initial conditions
     V0  = -70.
     Ca0 = 0.5
     x0  = [V0, mNa_inf(V0), hNa_inf(V0), mCaT_inf(V0), hCaT_inf(V0), mCaS_inf(V0), hCaS_inf(V0),
            mA_inf(V0), hA_inf(V0), mKCa_inf(V0, Ca0), mKd_inf(V0), mH_inf(V0), Ca0,
            gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak, gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak]

    # Simulation
    prob = ODEProblem(STG_homeo_leak_ODE, x0, (0, Tfinal/2), p) # Describing the problem
    sol = solve(prob) # Solving the problem

    return sol
end

## STG model from Liu 1998 - current-clamp mode
function simulate_STG_second_half(g_all, gCaS_init_crash, gA_init_crash, Iapp, tau_Na, tau_CaT, tau_CaS, tau_A, tau_KCa, tau_Kd, tau_H, tau_leak, Ca_tgt, C)
    ## Simulation of the model in current-clamp mode
    # Extracting the maximal ion channel conductances
    (gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak) = g_all

    # Parameter vector for simulations
    p = (Iapp, tau_Na, tau_CaT, tau_CaS, tau_A, tau_KCa, tau_Kd, tau_H, tau_leak, tau_g, Ca_tgt, C)

     # Initial conditions
     V0  = -70.
     Ca0 = 0.5
     x0  = [V0, mNa_inf(V0), hNa_inf(V0), mCaT_inf(V0), hCaT_inf(V0), mCaS_inf(V0), hCaS_inf(V0),
            mA_inf(V0), hA_inf(V0), mKCa_inf(V0, Ca0), mKd_inf(V0), mH_inf(V0), Ca0,
            gNa, gCaT, gCaS_init_crash, gA_init_crash, gKCa, gKd, gH, gleak, gNa, gCaT, gCaS_init_crash, gA_init_crash, gKCa, gKd, gH, gleak]

    # Simulation
    prob = ODEProblem(STG_homeo_leak_ODE, x0, (0, Tfinal/2), p) # Describing the problem
    sol = solve(prob) # Solving the problem

    return sol
end

## STG model from Liu 1998 - current-clamp mode
function simulate_STG_population(g_all_init, Iapp, tau_Na, tau_g, Ca_tgt,
                                 C, ICs_th_init, tt)
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
    tt_moving_average_ = window_size_s/2 : padding_s : Tfinal/2000 - window_size_s/2

    Ca_ma_matrix = zeros(ncells, 2*length(tt_moving_average_))

    @showprogress "Computing..." for i = 1 : ncells
        # Extracting the maximal ion channel conductances
        (gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak) = g_all_init[i, :]

        # Computing homeostasis time constants
        tau_CaT = tau_Na * gNa / gCaT
        tau_CaS = tau_Na * gNa / gCaS
        tau_A = tau_Na * gNa / gA
        tau_KCa = tau_Na * gNa / gKCa
        tau_Kd = tau_Na * gNa / gKd
        tau_H = tau_Na * gNa / gH
        tau_leak = tau_Na * gNa / gleak

        # Simulate to initial blockade conditions
        sol1 = simulate_STG_first_half(g_all_init[i, :], Iapp, tau_Na, tau_CaT, tau_CaS, tau_A, tau_KCa, tau_Kd, tau_H, tau_leak, Ca_tgt, C)

        # Retrieving variables
        x1 = sol1(tt[2:5000:Int(round(end/2))])

        # Pushing variables into matrices
        Ca_ma_matrix[i, 1:Int(round(end/2))] = moving_average(sol1(tt)[13, 1:Int(round(end/2))], window_size_, padding_)
        gNa_matrix[i, 1:Int(round(end/2))]   = x1[14, :]
        gCaT_matrix[i, 1:Int(round(end/2))]  = x1[15, :]
        gCaS_matrix[i, 1:Int(round(end/2))]  = x1[16, :]
        gA_matrix[i, 1:Int(round(end/2))]    = x1[17, :]
        gKCa_matrix[i, 1:Int(round(end/2))]  = x1[18, :]
        gKd_matrix[i, 1:Int(round(end/2))]   = x1[19, :]
        gH_matrix[i, 1:Int(round(end/2))]    = x1[20, :]

        (gCaS_init_crash, gA_init_crash) = DICs_gmax_neuromodCaSA(x1[14, end], x1[15, end],
            x1[19, end], x1[18, end], x1[20, end], x1[21, end], -8., 4., ICs_th_init[i, 1])

        # Simulate to initial blockade conditions
        g_all_init2 = [x1[14, end], x1[15, end], x1[16, end], x1[17, end], x1[18, end], x1[19, end], x1[20, end], x1[21, end]]
        sol2 = simulate_STG_second_half(g_all_init2, gCaS_init_crash, gA_init_crash, Iapp,
                                         tau_Na, tau_CaT, tau_CaS, tau_A, tau_KCa, tau_Kd, tau_H, tau_leak, Ca_tgt, C)

        # Retrieving variables
        x2 = sol2(tt[2:5000:Int(round(end/2))])

        # Pushing variables into matrices
        Ca_ma_matrix[i, Int(round(end/2))+1:end] = moving_average(sol2(tt)[13, 1:Int(round(end/2))], window_size_, padding_)
        gNa_matrix[i, Int(round(end/2))+1:end]   = x2[14, :]
        gCaT_matrix[i, Int(round(end/2))+1:end]  = x2[15, :]
        gCaS_matrix[i, Int(round(end/2))+1:end]  = x2[16, :]
        gA_matrix[i, Int(round(end/2))+1:end]    = x2[17, :]
        gKCa_matrix[i, Int(round(end/2))+1:end]  = x2[18, :]
        gKd_matrix[i, Int(round(end/2))+1:end]   = x2[19, :]
        gH_matrix[i, Int(round(end/2))+1:end]    = x2[20, :]
    end

    return gNa_matrix, gCaT_matrix, gCaS_matrix, gA_matrix, gKCa_matrix,
           gKd_matrix, gH_matrix, gleak_matrix, Ca_ma_matrix
end


## STG model from Liu 1998 - current-clamp mode
function simulate_STG_population_beautiful(g_all_init, Iapp, tau_Na, tau_g, Ca_tgt,
                                           C, ICs_th_init, tt)
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
    tt_moving_average_ = window_size_s/2 : padding_s : Tfinal/2000 - window_size_s/2

    Ca_ma_matrix = zeros(ncells, 2*length(tt_moving_average_))

    j_beautiful = 5

    @showprogress "Computing..." for i = 1 : ncells
        # Extracting the maximal ion channel conductances
        (gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak) = g_all_init[i, :]

        # Computing homeostasis time constants
        tau_CaT = tau_Na * gNa / gCaT
        tau_CaS = tau_Na * gNa / gCaS
        tau_A = tau_Na * gNa / gA
        tau_KCa = tau_Na * gNa / gKCa
        tau_Kd = tau_Na * gNa / gKd
        tau_H = tau_Na * gNa / gH
        tau_leak = tau_Na * gNa / gleak

        # Simulate to initial blockade conditions
        sol1 = simulate_STG_first_half(g_all_init[i, :], Iapp, tau_Na, tau_CaT, tau_CaS, tau_A, tau_KCa, tau_Kd, tau_H, tau_leak, Ca_tgt, C)

        # Retrieving variables
        x1 = sol1(tt[2:5000:Int(round(end/2))])

        # Pushing variables into matrices
        Ca_ma_matrix[i, 1:Int(round(end/2))] = moving_average(sol1(tt)[13, 1:Int(round(end/2))], window_size_, padding_)
        gNa_matrix[i, 1:Int(round(end/2))]   = x1[14, :]
        gCaT_matrix[i, 1:Int(round(end/2))]  = x1[15, :]
        gCaS_matrix[i, 1:Int(round(end/2))]  = x1[16, :]
        gA_matrix[i, 1:Int(round(end/2))]    = x1[17, :]
        gKCa_matrix[i, 1:Int(round(end/2))]  = x1[18, :]
        gKd_matrix[i, 1:Int(round(end/2))]   = x1[19, :]
        gH_matrix[i, 1:Int(round(end/2))]    = x1[20, :]
        gleak_matrix[i, 1:Int(round(end/2))] = x1[21, :]

        (gCaS_init_crash, gA_init_crash) = DICs_gmax_neuromodCaSA(gNa_matrix[i, j_beautiful], gCaT_matrix[i, j_beautiful],
            gKd_matrix[i, j_beautiful], gKCa_matrix[i, j_beautiful], gH_matrix[i, j_beautiful], gleak_matrix[i, j_beautiful], -8., 4., ICs_th_init[i, 1])

        # Simulate to initial blockade conditions
        g_all_init2 = [gNa_matrix[i, j_beautiful], gCaT_matrix[i, j_beautiful], gCaS_matrix[i, j_beautiful], gA_matrix[i, j_beautiful],
                       gKCa_matrix[i, j_beautiful], gKd_matrix[i, j_beautiful], gH_matrix[i, j_beautiful], gleak_matrix[i, j_beautiful]]
        sol2 = simulate_STG_second_half(g_all_init2, gCaS_init_crash, gA_init_crash, Iapp,
                                         tau_Na, tau_CaT, tau_CaS, tau_A, tau_KCa, tau_Kd, tau_H, tau_leak, Ca_tgt, C)

        # Retrieving variables
        x2 = sol2(tt[2:5000:Int(round(end/2))])

        # Pushing variables into matrices
        Ca_ma_matrix[i, Int(round(end/2))+1:end] = moving_average(sol2(tt)[13, 1:Int(round(end/2))], window_size_, padding_)
        gNa_matrix[i, Int(round(end/2))+1:end]   = x2[14, :] .* gNa_matrix[i, Int(round(end/2))] ./ gNa_matrix[i, j_beautiful]
        gCaT_matrix[i, Int(round(end/2))+1:end]  = x2[15, :] .* gCaT_matrix[i, Int(round(end/2))] ./ gCaT_matrix[i, j_beautiful]
        gCaS_matrix[i, Int(round(end/2))+1:end]  = x2[16, :] .* gCaS_matrix[i, Int(round(end/2))] ./ gCaS_matrix[i, j_beautiful]
        gA_matrix[i, Int(round(end/2))+1:end]    = x2[17, :] .* gA_matrix[i, Int(round(end/2))] ./ gA_matrix[i, j_beautiful]
        gKCa_matrix[i, Int(round(end/2))+1:end]  = x2[18, :] .* gKCa_matrix[i, Int(round(end/2))] ./ gKCa_matrix[i, j_beautiful]
        gKd_matrix[i, Int(round(end/2))+1:end]   = x2[19, :] .* gKd_matrix[i, Int(round(end/2))] ./ gKd_matrix[i, j_beautiful]
        gH_matrix[i, Int(round(end/2))+1:end]    = x2[20, :] .* gH_matrix[i, Int(round(end/2))] ./ gH_matrix[i, j_beautiful]
        gleak_matrix[i, Int(round(end/2))+1:end] = x2[21, :] .* gleak_matrix[i, Int(round(end/2))] ./ gleak_matrix[i, j_beautiful]
    end

    return gNa_matrix, gCaT_matrix, gCaS_matrix, gA_matrix, gKCa_matrix,
           gKd_matrix, gH_matrix, gleak_matrix, Ca_ma_matrix
end
