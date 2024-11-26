#=
This file contains differential equations describing the STG model
=#

include("network_STG_kinetics.jl") # Include STG model gating functions
using ProgressMeter

# STG ODEs
function dV(V, mNa, hNa, mCaT, hCaT, mCaS, hCaS, mA, hA, mKCa, mKd, mH,
            Ca, Iapp, gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak)
  (dt) * (1/C) * (- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                    gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) -
                    gKCa*mKCa^4*(V-VK) - gKd*mKd^4*(V-VK) - gH*mH*(V-VH) -
                    gleak*(V-Vleak) + Iapp)
end
dmNa(V, mNa) = (dt) * ((1/tau_mNa(V)) * (mNa_inf(V) - mNa))
dhNa(V, hNa) = (dt) * ((1/tau_hNa(V)) * (hNa_inf(V) - hNa))
dmCaT(V, mCaT) = (dt) * ((1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT))
dhCaT(V, hCaT) = (dt) * ((1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT))
dmCaS(V, mCaS) = (dt) * ((1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS))
dhCaS(V, hCaS) = (dt) * ((1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS))
dmA(V, mA) = (dt) * ((1/tau_mA(V)) * (mA_inf(V) - mA))
dhA(V, hA) = (dt) * ((1/tau_hA(V)) * (hA_inf(V) - hA))
dmKCa(V, Ca, mKCa, tmKCa) = (dt) * ((1/tau_mKCa(V, tmKCa)) * (mKCa_inf(V, Ca) - mKCa))
dmKd(V, mKd) = (dt) * ((1/tau_mKd(V)) * (mKd_inf(V) - mKd))
dmH(V, mH) = (dt) * ((1/tau_mH(V)) * (mH_inf(V) - mH))
dCa(V, mCaT, hCaT, mCaS, hCaS, Ca, gCaT, gCaS) = (dt) * (-0.94 * (gCaT*mCaT^3*hCaT*(V-VCa) + gCaS*mCaS^3*hCaS*(V-VCa)) - Ca + 0.05)/20
ds(V, s) = (dt) * ((1/tau_s(V)) * (s_inf(V) - s))
dgi(g, gi, α, β, u) = (dt) * (α * (g - gi) - β * gi + u)
dg(g, gi, α) = (dt) * (α * (gi - g))
dz(e, Kt, u, v) = (dt) * (e + Kt * (u - v))
dg_homeo(g, m, tau_g) = (dt) * ((1/tau_g) * (m - g))
dm_homeo(Ca, Ca_tgt, tau_X) = (dt) * ((1/tau_X) * (Ca_tgt - Ca))

# Control signal function
function uve(g, g_r, z, Kp, Ki, u_max)
    # Compute control signal
    e = g_r - g
    v = Kp * e + Ki * z

    # Anti-windup system
    if v > u_max
        u = u_max
    else
        u = v
    end

    return u, v, e
end

function simulateSTG_network(Iapp, gNainit, gCaTinit, gCaSinit, gAinit, gKCainit,
                             gKdinit, gHinit, gleakinit, tau_g, tau_Na, Ca_tgt,
                             gsyn21, gsyn12, gsyn13, gsyn53, gsyn54, gsyn45, gEl23, gEl43,
                             tmKCa, α, β, Kp, Ki, Kt, gsth, guth, Vth, u_maxCaS, u_maxA)

    # Initial conditions
    gNa = copy(gNainit)
    gCaT = copy(gCaTinit)
    gKCa = copy(gKCainit)
    gKd = copy(gKdinit)
    gH = copy(gHinit)
    gleak = copy(gleakinit)
    m_Na = copy(gNainit)
    m_CaT = copy(gCaTinit)
    m_KCa = copy(gKCainit)
    m_Kd = copy(gKdinit)
    m_H = copy(gHinit)
    m_leak = copy(gleakinit)
    tau_CaT = zero(tau_Na)
    tau_KCa = zero(tau_Na)
    tau_Kd = zero(tau_Na)
    tau_H = zero(tau_Na)
    tau_leak = zero(tau_Na)
    for j = 1 : 5
        tau_CaT[j] = tau_Na[j] * gNa[j] / gCaT[j]
        tau_KCa[j] = tau_Na[j] * gNa[j] / gKCa[j]
        tau_Kd[j] = tau_Na[j] * gNa[j] / gKd[j]
        tau_H[j] = tau_Na[j] * gNa[j] / gH[j]
        tau_leak[j] = tau_Na[j] * gNa[j] / gleak[j]
    end

    gCaS = copy(gCaSinit)
    gCaSi = copy(gCaSinit)
    zCaS = zero(gCaS)
    gA = copy(gAinit)
    gAi = copy(gAinit)
    zA = zero(gA)

    Vprev = .-60. .* ones(5) .+ 10. .* (rand(5).-0.5)
    Vcur = copy(Vprev)
    mNa = mNa_inf.(Vprev)
    hNa = hNa_inf.(Vprev)
    mCaT = mCaT_inf.(Vprev)
    hCaT = hCaT_inf.(Vprev)
    mCaS = mCaS_inf.(Vprev)
    hCaS = hCaS_inf.(Vprev)
    mA = mA_inf.(Vprev)
    hA = hA_inf.(Vprev)
    Ca = .-0.94 .* (gCaT .* mCaT.^3 .* hCaT .* (Vprev.-VCa) .+
                    gCaS .* mCaS.^3 .* hCaS .* (Vprev.-VCa)) .+ 0.05
    Caprev = copy(Ca)
    mKCa = mKCa_inf.(Vprev, Ca)
    mKd = mKd_inf.(Vprev)
    mH = mH_inf.(Vprev)

    gs_th = zeros(Tdt+1, 5)
    gu_th = guth[1](0.) # Assuming constant gu_th
    for j = 1 : 5
        gs_th[:, j] = gsth[j].(tsim)
    end

    gCaS_r = zeros(5)
    gA_r = zeros(5)
    for j = 1 : 5
        (gCaS_r[j], gA_r[j]) = DICs_gmax_neuromodCaSA_faster(gNa[j], gCaT[j], gKd[j], gKCa[j], gH[j],
                                                gleak[j], gs_th[1, j], gu_th, Vth[j], tmKCa[j])
    end


    s21 = s_inf(Vprev[2])
    s12 = s_inf(Vprev[1])
    s13 = s_inf(Vprev[1])
    s53 = s_inf(Vprev[5])
    s54 = s_inf(Vprev[5])
    s45 = s_inf(Vprev[4])

    # Initialize saving variables
    V_sol = zeros(length(tt_index), 5)
    Ca_sol = zero(V_sol)
    i = 1

    @showprogress "Computing ..." for z = 2 : Tdt+1
        # STG ODEs
        for j = 1 : 5
            # Liu model ODEs
            Vcur[j] += dV(Vprev[j], mNa[j], hNa[j], mCaT[j], hCaT[j], mCaS[j], hCaS[j], mA[j], hA[j], mKCa[j], mKd[j],
                          mH[j], Ca[j], Iapp[j], gNa[j], gCaT[j], gCaS[j], gA[j], gKCa[j], gKd[j], gH[j], gleak[j])
            mKCa[j] += dmKCa(Vprev[j], Ca[j], mKCa[j], tmKCa[j])
            Ca[j] += dCa(Vprev[j], mCaT[j], hCaT[j], mCaS[j], hCaS[j], Ca[j], gCaT[j], gCaS[j])
            mNa[j] += dmNa(Vprev[j], mNa[j])
            hNa[j] += dhNa(Vprev[j], hNa[j])
            mCaT[j] += dmCaT(Vprev[j], mCaT[j])
            hCaT[j] += dhCaT(Vprev[j], hCaT[j])
            mCaS[j] += dmCaS(Vprev[j], mCaS[j])
            hCaS[j] += dhCaS(Vprev[j], hCaS[j])
            mA[j] += dmA(Vprev[j], mA[j])
            hA[j] += dhA(Vprev[j], hA[j])
            mKd[j] += dmKd(Vprev[j], mKd[j])
            mH[j] += dmH(Vprev[j], mH[j])

            # Computing new reference values of gCaS and gA if DICs have changed
            (gCaS_r[j], gA_r[j]) = DICs_gmax_neuromodCaSA_faster(gNa[j], gCaT[j], gKd[j], gKCa[j], gH[j],
                                                gleak[j], gs_th[z, j], gu_th, Vth[j], tmKCa[j])

            # Computing control signals
            uCaS, vCaS, eCaS = uve(gCaS[j], gCaS_r[j], zCaS[j], Kp, Ki, u_maxCaS)
            uA, vA, eA = uve(gA[j], gA_r[j], zA[j], Kp, Ki, u_maxA)

            # Controller ODEs
            gCaSi[j] += dgi(gCaS[j], gCaSi[j], α, β, uCaS)
            gCaS[j] += dg(gCaS[j], gCaSi[j], α)
            zCaS[j] += dz(eCaS, Kt, uCaS, vCaS)
            gAi[j] += dgi(gA[j], gAi[j], α, β, uA)
            gA[j] += dg(gA[j], gAi[j], α)
            zA[j] += dz(eA, Kt, uA, vA)

            # Homeostasis ODEs
            gNa[j] += dg_homeo(gNa[j], m_Na[j], tau_g)
            gCaT[j] += dg_homeo(gCaT[j], m_CaT[j], tau_g)
            gKCa[j] += dg_homeo(gKCa[j], m_KCa[j], tau_g)
            gKd[j] += dg_homeo(gKd[j], m_Kd[j], tau_g)
            gH[j] += dg_homeo(gH[j], m_H[j], tau_g)
            gleak[j] += dg_homeo(gleak[j], m_leak[j], tau_g)

            m_Na[j] += dm_homeo(Caprev[j], Ca_tgt, tau_Na[j])
            m_CaT[j] += dm_homeo(Caprev[j], Ca_tgt, tau_CaT[j])
            m_KCa[j] += dm_homeo(Caprev[j], Ca_tgt, tau_KCa[j])
            m_Kd[j] += dm_homeo(Caprev[j], Ca_tgt, tau_Kd[j])
            m_H[j] += dm_homeo(Caprev[j], Ca_tgt, tau_H[j])
            m_leak[j] += dm_homeo(Caprev[j], Ca_tgt, tau_leak[j])

        end

        # Coupling ODEs
        Vcur[1] += (dt) * (1/C) * (-gsyn21*s21*(Vprev[1]-Vsyn))
        s21 += ds(Vprev[2], s21)

        Vcur[2] += (dt) * (1/C) * (-gsyn12*s12*(Vprev[2]-Vsyn) - gEl23*(Vprev[2]-Vprev[3]))
        s12 += ds(Vprev[1], s12)

        Vcur[3] += (dt) * (1/C) * (-gsyn13*s13*(Vprev[3]-Vsyn) - gsyn53*s53*(Vprev[3]-Vsyn) -
                                 gEl23*(Vprev[3]-Vprev[2]) - gEl43*(Vprev[3]-Vprev[4]))
        s13 += ds(Vprev[1], s13)
        s53 += ds(Vprev[5], s53)

        Vcur[4] += (dt) * (1/C) * (-gsyn54*s54*(Vprev[4]-Vsyn) - gEl43*(Vprev[4]-Vprev[3]))
        s54 += ds(Vprev[5], s54)

        Vcur[5] += (dt) * (1/C) * (-gsyn45*s45*(Vprev[5]-Vsyn))
        s45 += ds(Vprev[4], s45)

        Vprev = copy(Vcur)
        Caprev = copy(Ca)
        if (z-1) % dtratio == 0
            V_sol[i, :] = copy(Vcur)
            Ca_sol[i, :] = copy(Ca)
            i = i + 1
        end
    end

    return V_sol, Ca_sol
end

function simulateSTG_network_bad(Iapp, gNainit, gCaTinit, gCaSinit, gAinit, gKCainit,
                                 gKdinit, gHinit, gleakinit, tau_g, tau_Na, Ca_tgt,
                                 gsyn21, gsyn12, gsyn13, gsyn53, gsyn54, gsyn45, gEl23, gEl43,
                                 tmKCa, gsth, guth, Vth, tswitch)

    # Initial conditions
    gNa = copy(gNainit)
    gCaT = copy(gCaTinit)
    gCaS = copy(gCaSinit)
    gA = copy(gAinit)
    gKCa = copy(gKCainit)
    gKd = copy(gKdinit)
    gH = copy(gHinit)
    gleak = copy(gleakinit)
    m_Na = copy(gNainit)
    m_CaT = copy(gCaTinit)
    m_CaS = copy(gCaSinit)
    m_A = copy(gAinit)
    m_KCa = copy(gKCainit)
    m_Kd = copy(gKdinit)
    m_H = copy(gHinit)
    m_leak = copy(gleakinit)
    tau_CaT = zero(tau_Na)
    tau_CaS = zero(tau_Na)
    tau_A = zero(tau_Na)
    tau_KCa = zero(tau_Na)
    tau_Kd = zero(tau_Na)
    tau_H = zero(tau_Na)
    tau_leak = zero(tau_Na)
    for j = 1 : 5
        tau_CaT[j] = tau_Na[j] * gNa[j] / gCaT[j]
        tau_CaS[j] = tau_Na[j] * gNa[j] / gCaS[j]
        tau_A[j] = tau_Na[j] * gNa[j] / gA[j]
        tau_KCa[j] = tau_Na[j] * gNa[j] / gKCa[j]
        tau_Kd[j] = tau_Na[j] * gNa[j] / gKd[j]
        tau_H[j] = tau_Na[j] * gNa[j] / gH[j]
        tau_leak[j] = tau_Na[j] * gNa[j] / gleak[j]
    end

    Vprev = .-60. .* ones(5) .+ 10. .* (rand(5).-0.5)
    Vcur = copy(Vprev)
    mNa = mNa_inf.(Vprev)
    hNa = hNa_inf.(Vprev)
    mCaT = mCaT_inf.(Vprev)
    hCaT = hCaT_inf.(Vprev)
    mCaS = mCaS_inf.(Vprev)
    hCaS = hCaS_inf.(Vprev)
    mA = mA_inf.(Vprev)
    hA = hA_inf.(Vprev)
    Ca = .-0.94 .* (gCaT .* mCaT.^3 .* hCaT .* (Vprev.-VCa) .+
                    gCaS .* mCaS.^3 .* hCaS .* (Vprev.-VCa)) .+ 0.05
    Caprev = copy(Ca)
    mKCa = mKCa_inf.(Vprev, Ca)
    mKd = mKd_inf.(Vprev)
    mH = mH_inf.(Vprev)

    gs_th = gsth[5].(Tfinal)
    gu_th = guth[1](0.) # Assuming constant gu_th


    s21 = s_inf(Vprev[2])
    s12 = s_inf(Vprev[1])
    s13 = s_inf(Vprev[1])
    s53 = s_inf(Vprev[5])
    s54 = s_inf(Vprev[5])
    s45 = s_inf(Vprev[4])

    # Initialize saving variables
    V_sol = zeros(length(tt_index), 5)
    Ca_sol = zero(V_sol)
    i = 1

    @showprogress "Computing ..." for z = 2 : Tdt+1
        # Neuromodulating in the bad way
        if z == Int(tswitch/dt)
            for j = 4 : 5
                # Computing new reference values of gCaS and gA if DICs have changed
                (gCaS[j], gA[j]) = DICs_gmax_neuromodCaSA_faster(gNa[j], gCaT[j], gKd[j], gKCa[j], gH[j],
                                                    gleak[j], -8., 4., Vth[j], tmKCa[j])
            end
        end

        # STG ODEs
        for j = 1 : 5
            # Liu model ODEs
            Vcur[j] += dV(Vprev[j], mNa[j], hNa[j], mCaT[j], hCaT[j], mCaS[j], hCaS[j], mA[j], hA[j], mKCa[j], mKd[j],
                          mH[j], Ca[j], Iapp[j], gNa[j], gCaT[j], gCaS[j], gA[j], gKCa[j], gKd[j], gH[j], gleak[j])
            mKCa[j] += dmKCa(Vprev[j], Ca[j], mKCa[j], tmKCa[j])
            Ca[j] += dCa(Vprev[j], mCaT[j], hCaT[j], mCaS[j], hCaS[j], Ca[j], gCaT[j], gCaS[j])
            mNa[j] += dmNa(Vprev[j], mNa[j])
            hNa[j] += dhNa(Vprev[j], hNa[j])
            mCaT[j] += dmCaT(Vprev[j], mCaT[j])
            hCaT[j] += dhCaT(Vprev[j], hCaT[j])
            mCaS[j] += dmCaS(Vprev[j], mCaS[j])
            hCaS[j] += dhCaS(Vprev[j], hCaS[j])
            mA[j] += dmA(Vprev[j], mA[j])
            hA[j] += dhA(Vprev[j], hA[j])
            mKd[j] += dmKd(Vprev[j], mKd[j])
            mH[j] += dmH(Vprev[j], mH[j])

            # Homeostasis ODEs
            if (z > Int(2*tswitch/dt))
                gNa[j] += dg_homeo(gNa[j], m_Na[j], tau_g)
                gCaT[j] += dg_homeo(gCaT[j], m_CaT[j], tau_g)
                gCaS[j] += dg_homeo(gCaS[j], m_CaS[j], tau_g)
                gA[j] += dg_homeo(gA[j], m_A[j], tau_g)
                gKCa[j] += dg_homeo(gKCa[j], m_KCa[j], tau_g)
                gKd[j] += dg_homeo(gKd[j], m_Kd[j], tau_g)
                gH[j] += dg_homeo(gH[j], m_H[j], tau_g)
                gleak[j] += dg_homeo(gleak[j], m_leak[j], tau_g)

                m_Na[j] += dm_homeo(Caprev[j], Ca_tgt, tau_Na[j])
                m_CaT[j] += dm_homeo(Caprev[j], Ca_tgt, tau_CaT[j])
                m_CaS[j] += dm_homeo(Caprev[j], Ca_tgt, tau_CaS[j])
                m_A[j] += dm_homeo(Caprev[j], Ca_tgt, tau_A[j])
                m_KCa[j] += dm_homeo(Caprev[j], Ca_tgt, tau_KCa[j])
                m_Kd[j] += dm_homeo(Caprev[j], Ca_tgt, tau_Kd[j])
                m_H[j] += dm_homeo(Caprev[j], Ca_tgt, tau_H[j])
                m_leak[j] += dm_homeo(Caprev[j], Ca_tgt, tau_leak[j])
            end
        end

        # Coupling ODEs
        Vcur[1] += (dt) * (1/C) * (-gsyn21*s21*(Vprev[1]-Vsyn))
        s21 += ds(Vprev[2], s21)

        Vcur[2] += (dt) * (1/C) * (-gsyn12*s12*(Vprev[2]-Vsyn) - gEl23*(Vprev[2]-Vprev[3]))
        s12 += ds(Vprev[1], s12)

        Vcur[3] += (dt) * (1/C) * (-gsyn13*s13*(Vprev[3]-Vsyn) - gsyn53*s53*(Vprev[3]-Vsyn) -
                                 gEl23*(Vprev[3]-Vprev[2]) - gEl43*(Vprev[3]-Vprev[4]))
        s13 += ds(Vprev[1], s13)
        s53 += ds(Vprev[5], s53)

        Vcur[4] += (dt) * (1/C) * (-gsyn54*s54*(Vprev[4]-Vsyn) - gEl43*(Vprev[4]-Vprev[3]))
        s54 += ds(Vprev[5], s54)

        Vcur[5] += (dt) * (1/C) * (-gsyn45*s45*(Vprev[5]-Vsyn))
        s45 += ds(Vprev[4], s45)

        Vprev = copy(Vcur)
        Caprev = copy(Ca)
        if (z-1) % dtratio == 0
            V_sol[i, :] = copy(Vcur)
            Ca_sol[i, :] = copy(Ca)
            i = i + 1
        end
    end

    return V_sol, Ca_sol
end

## Stimulation function
heaviside(t) = (1 + sign(t)) / 2
pulse(t, ti, tf) = heaviside(t-ti) - heaviside(t-tf)
