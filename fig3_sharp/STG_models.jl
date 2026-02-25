#=
This file contains differential equations describing the STG model
=#

include("STG_kinetics.jl") # Include STG model gating functions

## STG model from Liu 1998 - current-clamp mode
function STG_homeo_leak_ODE(dx, x, p, t)
    # Parameters
    Iapp     = p[1]     # Amplitude of constant applied current
    tau_Na   = p[2]     # Sodium current time constant
    tau_CaT  = p[3]     # T-type calcium current time constant
    tau_CaS  = p[4]     # Slow calcium current time constant
    tau_A    = p[5]     # A-type potassium current time constant
    tau_KCa  = p[6]     # Calcium controlled potassium current time constant
    tau_Kd   = p[7]     # Delayed-rectifier potassium current time constant
    tau_H    = p[8]     # H current time constant
    tau_leak = p[9]     # Leak current time constant
    tau_g    = p[10]    # Translation time constant
    Ca_tgt   = p[11](t) # Calcium target
    C        = p[12]    # Membrane capacitance

    # Variables
    V      = x[1]  # Membrane potential
    mNa    = x[2]  # Sodium current activation
    hNa    = x[3]  # Sodium current inactivation
    mCaT   = x[4]  # T-type calcium current activation
    hCaT   = x[5]  # T-type calcium current inactivation
    mCaS   = x[6]  # Slow calcium current activation
    hCaS   = x[7]  # Slow calcium current inactivation
    mA     = x[8]  # A-type potassium current activation
    hA     = x[9]  # A-type potassium current inactivation
    mKCa   = x[10] # Calcium controlled potassium current activation
    mKd    = x[11] # Delayed-rectifier potassium current activation
    mH     = x[12] # H current activation
    Ca     = x[13] # Calcium concentration
    gNa    = x[14] # Sodium current maximal conductance
    gCaT   = x[15] # T-type calcium current maximal conductance
    gCaS   = x[16] # Slow calcium current maximal conductance
    gA     = x[17] # A-type potassium current maximal conductance
    gKCa   = x[18] # Calcium controlled potassium current maximal conductance
    gKd    = x[19] # Delayed-rectifier potassium current maximal conductance
    gH     = x[20] # H current maximal conductance
    gleak  = x[21] # Leak current maximal conductance
    m_Na   = x[22] # Sodium current mRNA
    m_CaT  = x[23] # T-type calcium current mRNA
    m_CaS  = x[24] # Slow calcium current mRNA
    m_A    = x[25] # A-type potassium current mRNA
    m_KCa  = x[26] # Calcium controlled potassium current mRNA
    m_Kd   = x[27] # Delayed-rectifier potassium current mRNA
    m_H    = x[28] # H current mRNA
    m_leak = x[29] # Leak current mRNA

    # ODEs
    dx[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - max(gA, 0)*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
                   gKd*mKd^4*(V-VK) - gH*mH*(V-VH) - gleak*(V-Vleak) + Iapp)

    dx[2] = (1/tau_mNa(V)) * (mNa_inf(V) - mNa)
    dx[3] = (1/tau_hNa(V)) * (hNa_inf(V) - hNa)
    dx[4] = (1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT)
    dx[5] = (1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT)
    dx[6] = (1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS)
    dx[7] = (1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS)
    dx[8] = (1/tau_mA(V)) * (mA_inf(V) - mA)
    dx[9] = (1/tau_hA(V)) * (hA_inf(V) - hA)
    dx[10] = (1/tau_mKCa(V)) * (mKCa_inf(V, Ca) - mKCa)
    dx[11] = (1/tau_mKd(V)) * (mKd_inf(V) - mKd)
    dx[12] = (1/tau_mH(V)) * (mH_inf(V) - mH)

    dx[13] = (-0.94*(gCaT*mCaT^3*hCaT*(V-VCa) + gCaS*mCaS^3*hCaS*(V-VCa)) - Ca + 0.05)/20

    # Homeostasis ODEs
    dx[14] = 1/tau_g * (m_Na - gNa)
    dx[15] = 1/tau_g * (m_CaT - gCaT)
    dx[16] = 1/tau_g * (m_CaS - gCaS)
    dx[17] = 1/tau_g * (m_A - gA)
    dx[18] = 1/tau_g * (m_KCa - gKCa)
    dx[19] = 1/tau_g * (m_Kd - gKd)
    dx[20] = 1/tau_g * (m_H - gH)
    dx[21] = 1/tau_g * (m_leak - gleak)

    dx[22] = 1/tau_Na * (Ca_tgt - Ca)
    dx[23] = 1/tau_CaT * (Ca_tgt - Ca)
    dx[24] = 1/tau_CaS * (Ca_tgt - Ca)
    dx[25] = 1/tau_A * (Ca_tgt - Ca)
    dx[26] = 1/tau_KCa * (Ca_tgt - Ca)
    dx[27] = 1/tau_Kd * (Ca_tgt - Ca)
    dx[28] = 1/tau_H * (Ca_tgt - Ca)
    dx[29] = 1/tau_leak * (Ca_tgt - Ca)
end

## Stimulation function
heaviside(t) = (1 + sign(t)) / 2
pulse(t, ti, tf) = heaviside(t-ti) - heaviside(t-tf)
