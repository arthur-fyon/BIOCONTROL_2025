### Kinetics and parameters for the DA model
# DA gating Functions
DA_boltz(V, A, B) = 1 / (1 + exp(-(V-A) / B))
tauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))

# Initializing Nernst reversalPotential
DA_VNa = 60. # Sodium reversal potential
DA_VK = -85. # Potassium reversal potential
DA_VCa = 60. # Calcium reversal potential
DA_VNMDA = 0. # NMDA reversal potential
DA_Vleak = -50. # Reversal potential of leak channels
DA_Mg = 1.4 # Mg concentration

# Na current
DA_mNa_inf(V) = DA_boltz(V, -30.0907, 9.7264)
DA_tau_mNa(V) = 0.01 + 1.0 / ((-(15.6504 + 0.4043*V)/(exp(-19.565 -0.5052*V)-1.0)) + 3.0212*exp(-7.4630e-3*V))
DA_hNa_inf(V) = DA_boltz(V, -54.0289, -10.7665)
DA_tau_hNa(V) = 0.4 + 1.0 / ((5.0754e-4*exp(-6.3213e-2*V)) + 9.7529*exp(0.13442*V))

# Kd current
DA_mKd_inf(V) = DA_boltz(V, -25., 12.)
DA_tau_mKd(V) = tauX(V, 20., 18., 38., -10.)

# CaL current
DA_mCaL_inf(V) = DA_boltz(V, -50., 2.)
DA_tau_mCaL(V) = tauX(V, 30., 28., 45., -3.)

# CaN current
DA_mCaN_inf(V) = DA_boltz(V, -30., 7.)
DA_tau_mCaN(V) = tauX(V, 30., 25., 55., -6.)

# ERG current
a0ERG(V) = 0.0036 * exp(0.0759*V)
b0ERG(V) = 1.2523e-5 * exp(-0.0671*V)
aiERG(V) = 0.1 * exp(0.1189*V)
biERG(V) = 0.003 * exp(-0.0733*V)
DA_o_inf(V) = a0ERG(V)*biERG(V) / (a0ERG(V)*(aiERG(V)+biERG(V)) + b0ERG(V)*biERG(V))
DA_tau_o(V) = 100000.

# NMDA current
DA_NMDA_inf(V, Mg) = 1 / (1 + Mg*exp(-0.08*V)/10.)
DA_tau_NMDA(V) = 1e-10
