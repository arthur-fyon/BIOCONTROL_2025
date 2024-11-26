#=
This file contains differential equations describing the STG model
=#

include("STG_kinetics.jl") # Include STG model gating functions

## STG model from Liu 1998 - current-clamp mode
function STG_ODE(du, u, p, t)
    # Parameters
    Iapp  = p[1]  # Amplitude of constant applied current
    gNa   = p[2]  # Sodium current maximal conductance
    gCaT  = p[3]  # T-type calcium current maximal conductance
    gCaS  = p[4]  # Slow calcium current maximal conductance
    gA    = p[5]  # A-type potassium current maximal conductance
    gKCa  = p[6]  # Calcium controlled potassium current maximal conductance
    gKd   = p[7]  # Delayed-rectifier potassium current maximal conductance
    gH    = p[8]  # H current maximal conductance
    gleak = p[9]  # Leak current maximal conductance
    C     = p[10] # Membrane capacitance

    # Variables
    V    = u[1]  # Membrane potential
    mNa  = u[2]  # Sodium current activation
    hNa  = u[3]  # Sodium current inactivation
    mCaT = u[4]  # T-type calcium current activation
    hCaT = u[5]  # T-type calcium current inactivation
    mCaS = u[6]  # Slow calcium current activation
    hCaS = u[7]  # Slow calcium current inactivation
    mA   = u[8]  # A-type potassium current activation
    hA   = u[9]  # A-type potassium current inactivation
    mKCa = u[10] # Calcium controlled potassium current activation
    mKd  = u[11] # Delayed-rectifier potassium current activation
    mH   = u[12] # H current activation
    Ca   = u[13] # Calcium concentration

    # ODEs
    du[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
                   gKd*mKd^4*(V-VK) - gH*mH*(V-VH) - gleak*(V-Vleak) + Iapp)

    du[2] = (1/tau_mNa(V)) * (mNa_inf(V) - mNa)
    du[3] = (1/tau_hNa(V)) * (hNa_inf(V) - hNa)
    du[4] = (1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT)
    du[5] = (1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT)
    du[6] = (1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS)
    du[7] = (1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS)
    du[8] = (1/tau_mA(V)) * (mA_inf(V) - mA)
    du[9] = (1/tau_hA(V)) * (hA_inf(V) - hA)
    du[10] = (1/tau_mKCa(V)) * (mKCa_inf(V, Ca) - mKCa)
    du[11] = (1/tau_mKd(V)) * (mKd_inf(V) - mKd)
    du[12] = (1/tau_mH(V)) * (mH_inf(V) - mH)

    du[13] = (-0.94*(gCaT*mCaT^3*hCaT*(V-VCa) + gCaS*mCaS^3*hCaS*(V-VCa)) - Ca + 0.05)/20
end

## STG model from Liu 1998 - current-clamp mode
function simple_PI_STG_ODE(dx, x, p, t)
    # Parameters
    Iapp     = p[1]  # Amplitude of constant applied current
    gNa      = p[2]  # Sodium current maximal conductance
    gCaT     = p[3]  # T-type calcium current maximal conductance
    gKCa     = p[4]  # Calcium controlled potassium current maximal conductance
    gKd      = p[5]  # Delayed-rectifier potassium current maximal conductance
    gH       = p[6]  # H current maximal conductance
    gleak    = p[7]  # Leak current maximal conductance
    C        = p[8]  # Membrane capacitance
    α        = p[9]  # Rate of transfer between intracellular and membrane
    β        = p[10] # Rate of degradation of intracellular proteins
    Kp       = p[11] # Proportional gain
    Ki       = p[12] # Integral gain
    Kt       = p[13] # Anti-windup gain
    gsth     = p[14] # Reference gs(Vth)
    guth     = p[15] # Reference gu(Vth)
    Vth      = p[16] # Threshold voltage
    u_maxCaS = p[17] # Maximal value of CaS actuator
    u_maxA   = p[18] # Maximal value of A actuator

    # Variables
    V        = x[1]     # Membrane potential
    mNa      = x[2]     # Sodium current activation
    hNa      = x[3]     # Sodium current inactivation
    mCaT     = x[4]     # T-type calcium current activation
    hCaT     = x[5]     # T-type calcium current inactivation
    mCaS     = x[6]     # Slow calcium current activation
    hCaS     = x[7]     # Slow calcium current inactivation
    mA       = x[8]     # A-type potassium current activation
    hA       = x[9]     # A-type potassium current inactivation
    mKCa     = x[10]    # Calcium controlled potassium current activation
    mKd      = x[11]    # Delayed-rectifier potassium current activation
    mH       = x[12]    # H current activation
    Ca       = x[13]    # Calcium concentration
    gCaSi    = x[14]    # Intracellular slow calcium current maximal conductance
    gCaS     = x[15]    # Slow calcium current maximal conductance
    zCaS     = x[16]    # Slow calcium current maximal conductance integral variable
    gAi      = x[17]    # Intracellular A-type potassium current maximal conductance
    gA       = x[18]    # A-type potassium current maximal conductance
    zA       = x[19]    # A-type potassium current maximal conductance integral variable

    # ODEs
    dx[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
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

    # Computing reference values of gCaS and gA
    gs_th = gsth(t)
    gu_th = guth(t)
    (gCaS_r, gA_r) = DICs_gmax_neuromodCaSA(gNa, gCaT, gKd, gKCa, gH, gleak,
                                              gs_th, gu_th, Vth)

    # Computing control signals
    eCaS = gCaS_r - gCaS
    vCaS = Kp * eCaS + Ki * zCaS
    eA = gA_r - gA
    vA = Kp * eA + Ki * zA

    # Anti-windup system
    if vCaS > u_maxCaS
        uCaS = u_maxCaS
    else
        uCaS = vCaS
    end
    if vA > u_maxA
        uA = u_maxA
    else
        uA = vA
    end

    # ODEs
    dx[14] = α * (gCaS - gCaSi) - β * gCaSi + uCaS
    dx[15] = α * (gCaSi - gCaS)
    dx[16] = eCaS + Kt * (uCaS - vCaS)
    dx[17] = α * (gA - gAi) - β * gAi + uA
    dx[18] = α * (gAi - gA)
    dx[19] = eA + Kt * (uA - vA)
end

## STG model from Liu 1998 - current-clamp mode
function STG_homeo_ODE(dx, x, p, t)
    # Parameters
    Iapp     = p[1]     # Amplitude of constant applied current
    tau_Na   = p[2]     # Sodium current time constant
    tau_CaT  = p[3]     # T-type calcium current time constant
    tau_CaS  = p[4]     # Slow calcium current time constant
    tau_A    = p[5]     # A-type potassium current time constant
    tau_KCa  = p[6]     # Calcium controlled potassium current time constant
    tau_Kd   = p[7]     # Delayed-rectifier potassium current time constant
    tau_H    = p[8]     # H current time constant
    gleak    = p[9]     # Leak current maximal conductance
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
    m_Na   = x[21] # Sodium current mRNA
    m_CaT  = x[22] # T-type calcium current mRNA
    m_CaS  = x[23] # Slow calcium current mRNA
    m_A    = x[24] # A-type potassium current mRNA
    m_KCa  = x[25] # Calcium controlled potassium current mRNA
    m_Kd   = x[26] # Delayed-rectifier potassium current mRNA
    m_H    = x[27] # H current mRNA

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

    dx[21] = 1/tau_Na * (Ca_tgt - Ca)
    dx[22] = 1/tau_CaT * (Ca_tgt - Ca)
    dx[23] = 1/tau_CaS * (Ca_tgt - Ca)
    dx[24] = 1/tau_A * (Ca_tgt - Ca)
    dx[25] = 1/tau_KCa * (Ca_tgt - Ca)
    dx[26] = 1/tau_Kd * (Ca_tgt - Ca)
    dx[27] = 1/tau_H * (Ca_tgt - Ca)
end

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

## STG model from Liu 1998 - current-clamp mode
function simple_PI_homeo_STG_ODE(dx, x, p, t)
    # Parameters
    Iapp     = p[1]     # Amplitude of constant applied current
    tau_Na   = p[2]     # Sodium current time constant
    tau_CaT  = p[3]     # T-type calcium current time constant
    tau_KCa  = p[4]     # Calcium controlled potassium current time constant
    tau_Kd   = p[5]     # Delayed-rectifier potassium current time constant
    tau_H    = p[6]     # H current time constant
    gleak    = p[7]     # Leak current maximal conductance
    tau_g    = p[8]     # Translation time constant
    Ca_tgt   = p[9](t)  # Calcium target
    C        = p[10]    # Membrane capacitance
    α        = p[11]    # Rate of transfer between intracellular and membrane
    β        = p[12]    # Rate of degradation of intracellular proteins
    Kp       = p[13]    # Proportional gain
    Ki       = p[14]    # Integral gain
    Kt       = p[15]    # Anti-windup gain
    gsth     = p[16](t) # Reference gs(Vth)
    guth     = p[17](t) # Reference gu(Vth)
    Vth      = p[18]    # Threshold voltage
    u_maxCaS = p[19]    # Maximal value of CaS actuator
    u_maxA   = p[20]    # Maximal value of A actuator

    # Variables
    V        = x[1]  # Membrane potential
    mNa      = x[2]  # Sodium current activation
    hNa      = x[3]  # Sodium current inactivation
    mCaT     = x[4]  # T-type calcium current activation
    hCaT     = x[5]  # T-type calcium current inactivation
    mCaS     = x[6]  # Slow calcium current activation
    hCaS     = x[7]  # Slow calcium current inactivation
    mA       = x[8]  # A-type potassium current activation
    hA       = x[9]  # A-type potassium current inactivation
    mKCa     = x[10] # Calcium controlled potassium current activation
    mKd      = x[11] # Delayed-rectifier potassium current activation
    mH       = x[12] # H current activation
    Ca       = x[13] # Calcium concentration
    gNa      = x[14] # Sodium current maximal conductance
    gCaT     = x[15] # T-type calcium current maximal conductance
    gKCa     = x[16] # Calcium controlled potassium current maximal conductance
    gKd      = x[17] # Delayed-rectifier potassium current maximal conductance
    gH       = x[18] # H current maximal conductance
    m_Na     = x[19] # Sodium current mRNA
    m_CaT    = x[20] # T-type calcium current mRNA
    m_KCa    = x[21] # Calcium controlled potassium current mRNA
    m_Kd     = x[22] # Delayed-rectifier potassium current mRNA
    m_H      = x[23] # H current mRNA
    gCaSi    = x[24] # Intracellular slow calcium current maximal conductance
    gCaS     = x[25] # Slow calcium current maximal conductance
    zCaS     = x[26] # Slow calcium current maximal conductance integral variable
    gAi      = x[27] # Intracellular A-type potassium current maximal conductance
    gA       = x[28] # A-type potassium current maximal conductance
    zA       = x[29] # A-type potassium current maximal conductance integral variable

    # ODEs
    dx[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
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
    dx[16] = 1/tau_g * (m_KCa - gKCa)
    dx[17] = 1/tau_g * (m_Kd - gKd)
    dx[18] = 1/tau_g * (m_H - gH)

    dx[19] = 1/tau_Na * (Ca_tgt - Ca)
    dx[20] = 1/tau_CaT * (Ca_tgt - Ca)
    dx[21] = 1/tau_KCa * (Ca_tgt - Ca)
    dx[22] = 1/tau_Kd * (Ca_tgt - Ca)
    dx[23] = 1/tau_H * (Ca_tgt - Ca)

    # Computing reference values of gCaS and gA
    (gCaS_r, gA_r) = DICs_gmax_neuromodCaSA(gNa, gCaT, gKd, gKCa, gH, gleak,
                                            gsth, guth, Vth)

    # Computing control signals
    eCaS = gCaS_r - gCaS
    vCaS = Kp * eCaS + Ki * zCaS
    eA = gA_r - gA
    vA = Kp * eA + Ki * zA

    # Anti-windup system
    if vCaS > u_maxCaS
        uCaS = u_maxCaS
    else
        uCaS = vCaS
    end
    if vA > u_maxA
        uA = u_maxA
    else
        uA = vA
    end

    # ODEs
    dx[24] = α * (gCaS - gCaSi) - β * gCaSi + uCaS
    dx[25] = α * (gCaSi - gCaS)
    dx[26] = eCaS + Kt * (uCaS - vCaS)
    dx[27] = α * (gA - gAi) - β * gAi + uA
    dx[28] = α * (gAi - gA)
    dx[29] = eA + Kt * (uA - vA)
end

## STG model from Liu 1998 - current-clamp mode
function simple_PI_homeo_leak_STG_ODE(dx, x, p, t)
    # Parameters
    Iapp     = p[1]     # Amplitude of constant applied current
    tau_Na   = p[2]     # Sodium current time constant
    tau_CaT  = p[3]     # T-type calcium current time constant
    tau_KCa  = p[4]     # Calcium controlled potassium current time constant
    tau_Kd   = p[5]     # Delayed-rectifier potassium current time constant
    tau_H    = p[6]     # H current time constant
    tau_leak = p[7]     # Leak current time constant
    tau_g    = p[8]     # Translation time constant
    Ca_tgt   = p[9](t)  # Calcium target
    C        = p[10]    # Membrane capacitance
    α        = p[11]    # Rate of transfer between intracellular and membrane
    β        = p[12]    # Rate of degradation of intracellular proteins
    Kp       = p[13]    # Proportional gain
    Ki       = p[14]    # Integral gain
    Kt       = p[15]    # Anti-windup gain
    gsth     = p[16](t) # Reference gs(Vth)
    guth     = p[17](t) # Reference gu(Vth)
    Vth      = p[18]    # Threshold voltage
    u_maxCaS = p[19]    # Maximal value of CaS actuator
    u_maxA   = p[20]    # Maximal value of A actuator

    # Variables
    V        = x[1]  # Membrane potential
    mNa      = x[2]  # Sodium current activation
    hNa      = x[3]  # Sodium current inactivation
    mCaT     = x[4]  # T-type calcium current activation
    hCaT     = x[5]  # T-type calcium current inactivation
    mCaS     = x[6]  # Slow calcium current activation
    hCaS     = x[7]  # Slow calcium current inactivation
    mA       = x[8]  # A-type potassium current activation
    hA       = x[9]  # A-type potassium current inactivation
    mKCa     = x[10] # Calcium controlled potassium current activation
    mKd      = x[11] # Delayed-rectifier potassium current activation
    mH       = x[12] # H current activation
    Ca       = x[13] # Calcium concentration
    gNa      = x[14] # Sodium current maximal conductance
    gCaT     = x[15] # T-type calcium current maximal conductance
    gKCa     = x[16] # Calcium controlled potassium current maximal conductance
    gKd      = x[17] # Delayed-rectifier potassium current maximal conductance
    gH       = x[18] # H current maximal conductance
    gleak    = x[19] # Leak current maximal conductance
    m_Na     = x[20] # Sodium current mRNA
    m_CaT    = x[21] # T-type calcium current mRNA
    m_KCa    = x[22] # Calcium controlled potassium current mRNA
    m_Kd     = x[23] # Delayed-rectifier potassium current mRNA
    m_H      = x[24] # H current mRNA
    m_leak   = x[25] # Leak current mRNA
    gCaSi    = x[26] # Intracellular slow calcium current maximal conductance
    gCaS     = x[27] # Slow calcium current maximal conductance
    zCaS     = x[28] # Slow calcium current maximal conductance integral variable
    gAi      = x[29] # Intracellular A-type potassium current maximal conductance
    gA       = x[30] # A-type potassium current maximal conductance
    zA       = x[31] # A-type potassium current maximal conductance integral variable

    # ODEs
    dx[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
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
    dx[16] = 1/tau_g * (m_KCa - gKCa)
    dx[17] = 1/tau_g * (m_Kd - gKd)
    dx[18] = 1/tau_g * (m_H - gH)
    dx[19] = 1/tau_g * (m_leak - gleak)

    dx[20] = 1/tau_Na * (Ca_tgt - Ca)
    dx[21] = 1/tau_CaT * (Ca_tgt - Ca)
    dx[22] = 1/tau_KCa * (Ca_tgt - Ca)
    dx[23] = 1/tau_Kd * (Ca_tgt - Ca)
    dx[24] = 1/tau_H * (Ca_tgt - Ca)
    dx[25] = 1/tau_leak * (Ca_tgt - Ca)

    # Computing reference values of gCaS and gA
    (gCaS_r, gA_r) = DICs_gmax_neuromodCaSA(gNa, gCaT, gKd, gKCa, gH, gleak,
                                            gsth, guth, Vth)

    # Computing control signals
    eCaS = gCaS_r - gCaS
    vCaS = Kp * eCaS + Ki * zCaS
    eA = gA_r - gA
    vA = Kp * eA + Ki * zA

    # Anti-windup system
    if vCaS > u_maxCaS
        uCaS = u_maxCaS
    else
        uCaS = vCaS
    end
    if vA > u_maxA
        uA = u_maxA
    else
        uA = vA
    end

    # ODEs
    dx[26] = α * (gCaS - gCaSi) - β * gCaSi + uCaS
    dx[27] = α * (gCaSi - gCaS)
    dx[28] = eCaS + Kt * (uCaS - vCaS)
    dx[29] = α * (gA - gAi) - β * gAi + uA
    dx[30] = α * (gAi - gA)
    dx[31] = eA + Kt * (uA - vA)
end

## Stimulation function
heaviside(t) = (1 + sign(t)) / 2
pulse(t, ti, tf) = heaviside(t-ti) - heaviside(t-tf)
