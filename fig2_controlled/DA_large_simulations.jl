#=
This file contains large simulations scripts
=#

using DifferentialEquations, ProgressMeter

include("DA_kinetics.jl") # Include STG kinetics of gating variables
include("DA_ODE.jl") # Include STG model
include("STG_utils.jl") # Include some utils functions

# Moving average function
moving_average(vs, n, padding) = [sum(vs[i:(i+n-1)])/n for i in 1:padding:(length(vs)-(n-1))];

## STG model from Liu 1998 - current-clamp mode
function simulate_DA(g_all, Iapp, tau_Na, tau_g, Ca_tgt, C, Mg, α, β, Kp, Ki,
                     SVth, gsth_sim, guth_sim)
    ## Simulation of the model in current-clamp mode
    # Extracting the maximal ion channel conductances
    (gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak) = g_all

    # Computing homeostasis time constants
    tau_Kd = tau_Na * gNa / gKd
    tau_ERG = tau_Na * gNa / gERG
    tau_NMDA = tau_Na * gNa / gNMDA
    tau_leak = tau_Na * gNa / gleak

    # Parameter vector for simulations
    p = [Iapp, tau_Na, tau_Kd, tau_ERG, tau_NMDA, tau_leak, tau_g, Ca_tgt, C, Mg, α, β,
         Kp, Ki, SVth, gsth_sim, guth_sim]

     # Initial conditions
     V0 = -90.
     x0  = [V0, DA_mNa_inf(V0), DA_hNa_inf(V0), DA_mKd_inf(V0), DA_mCaL_inf(V0), DA_mCaN_inf(V0), 0., 0.,
            0.5, gNa, gKd, gERG, gNMDA, gleak, gNa, gKd, gERG, gNMDA, gleak, gCaL, gCaL, (β * gCaL) / Ki,
            gCaN, gCaN, (β * gCaN) / Ki]

    # Simulation
    prob = ODEProblem(DA_ODE_PI_homeo, x0, tspan, p) # Describing the problem
    sol = solve(prob) # Solving the problem

    return sol
end

## STG model from Liu 1998 - current-clamp mode
function simulate_DA_population(g_all_init, Iapp, tau_Na, tau_g, Ca_tgt,
                                 C, Mg, α, β, Kp, Ki, SVth,
                                 gsth_sim, guth_sim, tt)
    # Initializing some variables
    gNa_matrix = zeros(ncells, length(tt[2:5000:end]))
    gKd_matrix = zero(gNa_matrix)
    gCaL_matrix = zero(gNa_matrix)
    gCaN_matrix = zero(gNa_matrix)
    gERG_matrix = zero(gNa_matrix)
    gNMDA_matrix = zero(gNa_matrix)
    gleak_matrix = zero(gNa_matrix)

    window_size_s = 5
    padding_s = 0.2
    window_size_ = Int64(round(window_size_s*1000*length(tt)/(Tfinal)))
    padding_ = Int64(round(padding_s*1000*length(tt)/(Tfinal)))
    tt_moving_average_ = window_size_s/2 : padding_s : Tfinal/1000 - window_size_s/2

    Ca_ma_matrix = zeros(ncells, length(tt_moving_average_))

    @showprogress "Computing..." for i = 1 : ncells
        # Simulate to initial blockade conditions
        sol = simulate_DA(g_all_init[i, :], Iapp, tau_Na, tau_g, Ca_tgt, C, Mg, α, β, Kp, Ki,
                          SVth, gsth_sim, guth_sim)

        # Retrieving variables
        x = sol(tt[2:5000:end])

        # Pushing variables into matrices
        gNa_matrix[i, :] = x[10, :]
        gKd_matrix[i, :] = x[11, :]
        gCaL_matrix[i, :] = x[21, :]
        gCaN_matrix[i, :] = x[24, :]
        gERG_matrix[i, :] = x[12, :]
        gNMDA_matrix[i, :] = x[13, :]
        gleak_matrix[i, :] = x[14, :]

        Ca_ma_matrix[i, :] = moving_average(sol(tt)[9, :], window_size_, padding_)
    end

    return gNa_matrix, gKd_matrix, gCaL_matrix, gCaN_matrix,
           gERG_matrix, gNMDA_matrix, gleak_matrix, Ca_ma_matrix
end
