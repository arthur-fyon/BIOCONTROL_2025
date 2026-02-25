#=
This file contains differential equations describing the STG model
=#

include("STG_kinetics.jl") # Include STG model gating functions

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

## STG model from Liu 1998 - current-clamp mode
function CaT_blockade_STG_ODE(dx, x, p, t)
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
    gCaT     = 0        # Blockade of CaT currents

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
    gKCa     = x[15] # Calcium controlled potassium current maximal conductance
    gKd      = x[16] # Delayed-rectifier potassium current maximal conductance
    gH       = x[17] # H current maximal conductance
    gleak    = x[18] # Leak current maximal conductance
    m_Na     = x[19] # Sodium current mRNA
    m_KCa    = x[20] # Calcium controlled potassium current mRNA
    m_Kd     = x[21] # Delayed-rectifier potassium current mRNA
    m_H      = x[22] # H current mRNA
    m_leak   = x[23] # Leak current mRNA
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
    dx[15] = 1/tau_g * (m_KCa - gKCa)
    dx[16] = 1/tau_g * (m_Kd - gKd)
    dx[17] = 1/tau_g * (m_H - gH)
    dx[18] = 1/tau_g * (m_leak - gleak)

    dx[19] = 1/tau_Na * (Ca_tgt - Ca)
    dx[20] = 1/tau_KCa * (Ca_tgt - Ca)
    dx[21] = 1/tau_Kd * (Ca_tgt - Ca)
    dx[22] = 1/tau_H * (Ca_tgt - Ca)
    dx[23] = 1/tau_leak * (Ca_tgt - Ca)

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
function H_blockade_STG_ODE(dx, x, p, t)
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
    gH       = 0        # Blockade of H currents

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
    gleak    = x[18] # Leak current maximal conductance
    m_Na     = x[19] # Sodium current mRNA
    m_CaT    = x[20] # T-type calcium current mRNA
    m_KCa    = x[21] # Calcium controlled potassium current mRNA
    m_Kd     = x[22] # Delayed-rectifier potassium current mRNA
    m_leak   = x[23] # Leak current mRNA
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
    dx[18] = 1/tau_g * (m_leak - gleak)

    dx[19] = 1/tau_Na * (Ca_tgt - Ca)
    dx[20] = 1/tau_CaT * (Ca_tgt - Ca)
    dx[21] = 1/tau_KCa * (Ca_tgt - Ca)
    dx[22] = 1/tau_Kd * (Ca_tgt - Ca)
    dx[23] = 1/tau_leak * (Ca_tgt - Ca)

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
function Na_blockade_STG_ODE(dx, x, p, t)
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
    gNa      = 0        # Blockade of ion channels

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
    gCaT     = x[14] # T-type calcium current maximal conductance
    gKCa     = x[15] # Calcium controlled potassium current maximal conductance
    gKd      = x[16] # Delayed-rectifier potassium current maximal conductance
    gH       = x[17] # H current maximal conductance
    gleak    = x[18] # Leak current maximal conductance
    m_CaT    = x[19] # T-type calcium current mRNA
    m_KCa    = x[20] # Calcium controlled potassium current mRNA
    m_Kd     = x[21] # Delayed-rectifier potassium current mRNA
    m_H      = x[22] # H current mRNA
    m_leak   = x[23] # Leak current mRNA
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
    dx[14] = 1/tau_g * (m_CaT - gCaT)
    dx[15] = 1/tau_g * (m_KCa - gKCa)
    dx[16] = 1/tau_g * (m_Kd - gKd)
    dx[17] = 1/tau_g * (m_H - gH)
    dx[18] = 1/tau_g * (m_leak - gleak)

    dx[19] = 1/tau_CaT * (Ca_tgt - Ca)
    dx[20] = 1/tau_KCa * (Ca_tgt - Ca)
    dx[21] = 1/tau_Kd * (Ca_tgt - Ca)
    dx[22] = 1/tau_H * (Ca_tgt - Ca)
    dx[23] = 1/tau_leak * (Ca_tgt - Ca)

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
function KCa_blockade_STG_ODE(dx, x, p, t)
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
    gKCa     = 0        # Blockade of Ca-controlled potassium current

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
    gKd      = x[16] # Delayed-rectifier potassium current maximal conductance
    gH       = x[17] # H current maximal conductance
    gleak    = x[18] # Leak current maximal conductance
    m_Na     = x[19] # Sodium current mRNA
    m_CaT    = x[20] # T-type calcium current mRNA
    m_Kd     = x[21] # Delayed-rectifier potassium current mRNA
    m_H      = x[22] # H current mRNA
    m_leak   = x[23] # Leak current mRNA
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
    dx[16] = 1/tau_g * (m_Kd - gKd)
    dx[17] = 1/tau_g * (m_H - gH)
    dx[18] = 1/tau_g * (m_leak - gleak)

    dx[19] = 1/tau_Na * (Ca_tgt - Ca)
    dx[20] = 1/tau_CaT * (Ca_tgt - Ca)
    dx[21] = 1/tau_Kd * (Ca_tgt - Ca)
    dx[22] = 1/tau_H * (Ca_tgt - Ca)
    dx[23] = 1/tau_leak * (Ca_tgt - Ca)

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

## Stimulation function
heaviside(t) = (1 + sign(t)) / 2
pulse(t, ti, tf) = heaviside(t-ti) - heaviside(t-tf)
