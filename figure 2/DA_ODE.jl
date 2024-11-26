#=
	This file contains differential equations describing the controlled CB model of interest
=#

# Function that outputs values of variables derivatives
function DA_ODE_PI_homeo(dx, x, p, t)
	# Parameters
	Iapp 	 = p[1](t)     # Time dependent applied current
	tau_Na   = p[2]        # Sodium current time constant
    tau_Kd   = p[3]        # Delayed-rectifier potassium current time constant
    tau_ERG  = p[4]        # ERG current time constant
    tau_NMDA = p[5]        # NMDA current time constant
    tau_leak = p[6]        # Leak current time constant
    tau_g    = p[7]        # Translation time constant
    Ca_tgt   = p[8](t)     # Calcium target
	C 		 = p[9] 	   # Membrane capacitance
	Mg 		 = p[10]	   # Mg concentration
	α 		 = p[11] 	   # Rate of transfer between intracellular and membrane
	β 		 = p[12] 	   # Rate of degradation of intracellular proteins
	Kp 		 = p[13] 	   # Proportional gain
	Ki 		 = p[14] 	   # Integral gain
	SVth 	 = p[15]/x[14] # Sensitivity matrix at threshold voltage normalized
	gsth 	 = p[16](t)    # Reference gs(Vth)
	guth 	 = p[17](t)    # Reference gu(Vth)

	# Variables
	V      = x[1]  # Membrane voltage
	mNa    = x[2]  # Activation gating variable of current Na
	hNa    = x[3]  # Inactivation gating variable of current Na
	mKd    = x[4]  # Activation gating variable of current Kd
	mCaL   = x[5]  # Activation gating variable of current CaL
	mCaN   = x[6]  # Activation gating variable of current CaN
	oERG   = x[7]  # ERG potassium current activation
	iERG   = x[8]  # ERG current intermediate
	Ca     = x[9]  # Intracellular calcium
	gNa    = x[10] # Maximum conductance of current Na
	gKd    = x[11] # Maximum conductance of current Kd
	gERG   = x[12] # Maximum conductance of current ERG
	gNMDA  = x[13] # Maximum conductance of current NMDA
	gleak  = x[14] # Leakage conductance
	m_Na   = x[15] # Maximum conductance of current Na mRNA
	m_Kd   = x[16] # Maximum conductance of current Kd mRNA
	m_ERG  = x[17] # Maximum conductance of current ERG mRNA
	m_NMDA = x[18] # Maximum conductance of current NMDA mRNA
	m_leak = x[19] # Leakage conductance mRNA
	gCaLi  = x[20] # Intracellular maximum conductance of current CaL
	gCaL   = x[21] # Maximum conductance of current CaL
	zCaL   = x[22] # Integral variable of current CaL
	gCaNi  = x[23] # Intracellular maximum conductance of current CaN
	gCaN   = x[24] # Maximum conductance of current CaN
	zCaN   = x[25] # Integral variable of current CaN


	# ODEs
	dx[1] = (1/C) * (- gNa*mNa^3*hNa^1*(V - 60.0) -
					   gKd*mKd^3*(V - -85.0) -
					   gCaL*mCaL^2*(V - 60.0) -
					   gCaN*mCaN^1*(V - 60.0) -
					   gERG*oERG^1*(V - -85.0) -
					   gNMDA*DA_NMDA_inf(V, Mg)^1*(V - 0.0) -
					   gleak*(V - -50.0) + Iapp)
	dx[2] = (1/DA_tau_mNa(V)) * (DA_mNa_inf(V) - mNa)
	dx[3] = (1/DA_tau_hNa(V)) * (DA_hNa_inf(V) - hNa)
	dx[4] = (1/DA_tau_mKd(V)) * (DA_mKd_inf(V) - mKd)
	dx[5] = (1/DA_tau_mCaL(V)) * (DA_mCaL_inf(V) - mCaL)
	dx[6] = (1/DA_tau_mCaN(V)) * (DA_mCaN_inf(V) - mCaN)
	dx[7] = a0ERG(V) * (1-oERG-iERG) + biERG(V) * iERG - oERG * (aiERG(V)+b0ERG(V))
    dx[8] = aiERG(V) * oERG - biERG(V) * iERG

	# Added calcium dynamics
	dx[9] = (-0.94*gCaL*mCaL^2*(V - 60.0) + -0.94*gCaN*mCaN^1*(V - 60.0) + -Ca + 0.05) / 20.0

	# Homeostasis ODEs
    dx[10] = 1/tau_g * (m_Na - gNa)
    dx[11] = 1/tau_g * (m_Kd - gKd)
    dx[12] = 1/tau_g * (m_ERG - gERG)
    dx[13] = 1/tau_g * (m_NMDA - gNMDA)
    dx[14] = 1/tau_g * (m_leak - gleak)

    dx[15] = 1/tau_Na * (Ca_tgt - Ca)
    dx[16] = 1/tau_Kd * (Ca_tgt - Ca)
    dx[17] = 1/tau_ERG * (Ca_tgt - Ca)
    dx[18] = 1/tau_NMDA * (Ca_tgt - Ca)
    dx[19] = 1/tau_leak * (Ca_tgt - Ca)

	# Retrieving which line of the sensitivity matrix matter
	timescales = [2, 3]

	# Retrieving which column of the sensitivity matrix belong to unmodulated conductances
	unmodulated = [1, 2, 5, 6]

	# Retrieving which column of the sensitivity matrix belong to modulated conductances
	modulated = [3, 4]

	# Computing the right hand side of the linear system
	gDICr = [gsth, guth]
	gDICr = gDICr - SVth[timescales, unmodulated] * collect(x[10:13])

	# Computing the left hand side of the linear system
	Smod = SVth[timescales, modulated]

	# Computing the solution of the linear system
	g_r = \(Smod, gDICr)

	# Error signals and control inputs
	eCaL = g_r[1] - gCaL
	uCaL = Kp * eCaL + Ki * zCaL
	eCaN = g_r[2] - gCaN
	uCaN = Kp * eCaN + Ki * zCaN

	# ODEs of the controller
	dx[20] = α * gCaL - α * gCaLi - β * gCaLi + uCaL
	dx[21] = α * gCaLi - α * gCaL
	dx[22] = eCaL
	dx[23] = α * gCaN - α * gCaNi - β * gCaNi + uCaN
	dx[24] = α * gCaNi - α * gCaN
	dx[25] = eCaN

end
