#=
This file animates the simple PI on 1 neuron
=#

using Plots, LaTeXStrings, Random, GLMakie
using DataStructures: CircularBuffer
include("network_STG_kinetics.jl") # Loading of STG kinetics of gating variables
include("network_STG_models.jl") # Loading of STG model
include("network_STG_utils.jl") # Loading of some utils functions
include("network_STG_gs_derivatives.jl") # Loading of X_inf derivatives
include("network_STG_DIC.jl") # Loading of the DIC and compensation algorithm
include("network_STG_neuromodulation.jl"); # Loading of the neuromodulation cells functions

# Definition of simulation time (in ms)
const Tfinal = 20000
const dt = 0.0005
const tsim = 0 : dt : Tfinal
const dtplot = 0.2
const tt = 0 : dtplot : Tfinal
const dtratio = Int(dtplot/dt)
const Tdt = Int(Tfinal/dt)
const tt_index = 1 : dtratio : Tdt+1

# Definition of reversal potential values (in mV) and membrane capacitance
const VNa   = 50. # Sodium reversal potential
const VK    = -80. # Potassium reversal potential
const VCa   = 80. # Calcium reversal potential
const VH    = -20. # Reversal potential for the H-current (permeable to both sodium and potassium ions)
const Vleak = -50. # Reversal potential of leak channels
const Vsyn = -75. # Reversal potential of synaptic channels
const C     = 1. # Membrane capacitance

# Definition of voltage range for the DICs
const Vmin = -60
const Vmax = 0
const V    = range(Vmin, stop=Vmax, step=0.01)

# Definition of the number of cells in the random set
const ncells = 5

# Modifying backend GR attributes
gr(guidefontsize=18, legendfontsize=12, margin=5Plots.mm, grid=false)
default(fmt = :png)

# Initial firing pattern
guth = 4.
tmKCa = [2.; 2.; 2.; 20.; 20.]
(g_all_init, ICs_th_init) = degeneracy_fixDICs_neuromod(ncells, 5., guth, -50., tmKCa)
# create a spiking set with max variability in gCaS and gA

# Extracting the maximal ion channel conductances
gNa = g_all_init[:, 1]
gCaT = g_all_init[:, 2]
gCaSinit = g_all_init[:, 3]
gAinit = g_all_init[:, 4]
gKCa = g_all_init[:, 5]
gKd = g_all_init[:, 6]
gH = g_all_init[:, 7]
gleak = g_all_init[:, 8]

# Definition of parameters
α = 5e-3 # Rate of transfer between intracellular and membrane
β = 5e-3 # Rate of degradation of intracellular proteins
Kp = 3e-4 # Proprtional gain
Ki = 5e-6 # Integral gain
Kt = β / Ki # Anti-windup gain
gsth_sim = -8.
gs_th = [gsth_sim, gsth_sim, gsth_sim, gsth_sim, gsth_sim]
guth_sim = 4.
gu_th = [guth_sim, guth_sim, guth_sim, guth_sim, guth_sim]
u_max = 1e7 # Maximum value of actuator
Vth = ICs_th_init[:, 1]

gsyn21 = 0.2*4
gsyn12 = 0.2*4
gsyn13 = 0.2*4
gsyn53 = 0.2*4
gsyn54 = 0.2*4
gsyn45 = 0.2*4
gEl23 = 0.05
gEl43 = 0.05

# Input current definition
Iapp = 0. * ones(5)

# Function to initialize the animated figure
function makefig(Vcur, t)
    # Initializing the observable buffer containing the trace to be plotted
    tail = Int64(5000/0.2)
    trace1 = CircularBuffer{Point2f}(tail)
    trace2 = CircularBuffer{Point2f}(tail)
    trace3 = CircularBuffer{Point2f}(tail)
    trace4 = CircularBuffer{Point2f}(tail)
    trace5 = CircularBuffer{Point2f}(tail)
    fill!(trace1, Point2f(t[1], Vcur[1]))
    fill!(trace2, Point2f(t[1], Vcur[2]))
    fill!(trace3, Point2f(t[1], Vcur[3]))
    fill!(trace4, Point2f(t[1], Vcur[4]))
    fill!(trace5, Point2f(t[1], Vcur[5]))

    trace1 = Observable(trace1)
    trace2 = Observable(trace2)
    trace3 = Observable(trace3)
    trace4 = Observable(trace4)
    trace5 = Observable(trace5)

    # Creating the figure
    fig = GLMakie.Figure(resolution=(2000, 1500))

    # Customizing the axes
    ax1 = GLMakie.Axis(fig[1, 1:30], title="Neuron recordings", xlabel="", ylabel="",
                       xticklabelsize=30, yticklabelsize=30, titlesize=60, xlabelsize=50,
                       ylabelsize=50, xgridwidth=0, ygridwidth=0, yticks=[-100., 60.])
    GLMakie.xlims!(ax1, 0., 5.)
    GLMakie.ylims!(ax1, -100., 60.)

    ax2 = GLMakie.Axis(fig[2, 1:30], title="", xlabel="", ylabel="",
                       xticklabelsize=30, yticklabelsize=30, titlesize=60, xlabelsize=50,
                       ylabelsize=50, xgridwidth=0, ygridwidth=0, yticks=[-100., 60.])
    GLMakie.xlims!(ax2, 0., 5.)
    GLMakie.ylims!(ax2, -100., 60.)

    ax3 = GLMakie.Axis(fig[3, 1:30], title="", xlabel="", ylabel="V (mV)",
                       xticklabelsize=30, yticklabelsize=30, titlesize=60, xlabelsize=50,
                       ylabelsize=50, xgridwidth=0, ygridwidth=0, yticks=[-100., 60.])
    GLMakie.xlims!(ax3, 0., 5.)
    GLMakie.ylims!(ax3, -100., 60.)

    ax4 = GLMakie.Axis(fig[4, 1:30], title="", xlabel="", ylabel="",
                       xticklabelsize=30, yticklabelsize=30, titlesize=60, xlabelsize=50,
                       ylabelsize=50, xgridwidth=0, ygridwidth=0, yticks=[-100., 60.])
    GLMakie.xlims!(ax4, 0., 5.)
    GLMakie.ylims!(ax4, -100., 60.)

    ax5 = GLMakie.Axis(fig[5, 1:30], title="", xlabel="", ylabel="",
                       xticklabelsize=30, yticklabelsize=30, titlesize=60, xlabelsize=50,
                       ylabelsize=50, xgridwidth=0, ygridwidth=0, yticks=[-100., 60.])
    GLMakie.xlims!(ax5, 0., 5.)
    GLMakie.ylims!(ax5, -100., 60.)

    # Plotting the empty trace with the initial color of Julia with fading
    logocolors = Colors.JULIA_LOGO_COLORS
    c = to_color(logocolors.blue)
    tailcol = [RGBAf(c.r, c.g, c.b, (i/tail)^2) for i in 1:tail]
    GLMakie.lines!(ax1, trace1; linewidth = 3, color = tailcol)
    GLMakie.lines!(ax2, trace2; linewidth = 3, color = tailcol)
    GLMakie.lines!(ax3, trace3; linewidth = 3, color = tailcol)
    GLMakie.lines!(ax4, trace4; linewidth = 3, color = tailcol)
    GLMakie.lines!(ax5, trace5; linewidth = 3, color = tailcol)

    # Resizing the window
    screen = display(fig)
    GLMakie.resize!(screen, 2000, 1500)

    return tail, trace1, trace2, trace3, trace4, trace5, fig, screen
end

function progress_for_one_step!(Vcur, mNa, hNa, mCaT, hCaT, mCaS, hCaS, mA, hA,
                                mKCa, mKd, mH, Ca, gCaS, gCaSi, zCaS, gA, gAi, zA, gCaS_r, gA_r,
                                Iapp, gNa, gCaT, gKCa, gKd, gH, gleak, gsyn21, gsyn12,
                                gsyn13, gsyn53, gsyn54, gsyn45, gEl23, gEl43, tmKCa,
                                sxx, C, α, β, Kp, Ki, Kt, gsth, gsth_prev, guth, Vth, u_max, t)
    # Initializing variables from last time step
    Vprev = copy(Vcur)

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

        # Computing new reference values of gCaS and gA
        if !(gsth[j] == gsth_prev[j])
            (gCaS_r[j], gA_r[j]) = DICs_gmax_neuromodCaSA(gNa[j], gCaT[j], gKd[j], gKCa[j], gH[j],
                                            gleak[j], gsth[j], guth[j], Vth[j], tmKCa[j])
            gsth_prev[j] = gsth[j]
        end

        # Computing control signals
        uCaS, vCaS, eCaS = uve(gCaS[j], gCaS_r[j], zCaS[j], Kp, Ki, u_max)
        uA, vA, eA = uve(gA[j], gA_r[j], zA[j], Kp, Ki, u_max)

        # Controller ODEs
        gCaSi[j] += dgi(gCaS[j], gCaSi[j], α, β, uCaS)
        gCaS[j] += dg(gCaS[j], gCaSi[j], α)
        zCaS[j] += dz(eCaS, Kt, uCaS, vCaS)
        gAi[j] += dgi(gA[j], gAi[j], α, β, uA)
        gA[j] += dg(gA[j], gAi[j], α)
        zA[j] += dz(eA, Kt, uA, vA)
    end

    # Coupling ODEs
    Vcur[1] += (dt) * (1/C) * (-gsyn21*sxx[1]*(Vprev[1]-Vsyn))
    sxx[1] += ds(Vprev[2], sxx[1])

    Vcur[2] += (dt) * (1/C) * (-gsyn12*sxx[2]*(Vprev[2]-Vsyn) - gEl23*(Vprev[2]-Vprev[3]))
    sxx[2] += ds(Vprev[1], sxx[2])

    Vcur[3] += (dt) * (1/C) * (-gsyn13*sxx[3]*(Vprev[3]-Vsyn) - gsyn53*sxx[4]*(Vprev[3]-Vsyn) -
                             gEl23*(Vprev[3]-Vprev[2]) - gEl43*(Vprev[3]-Vprev[4]))
    sxx[3] += ds(Vprev[1], sxx[3])
    sxx[4] += ds(Vprev[5], sxx[4])

    Vcur[4] += (dt) * (1/C) * (-gsyn54*sxx[5]*(Vprev[4]-Vsyn) - gEl43*(Vprev[4]-Vprev[3]))
    sxx[5] += ds(Vprev[5], sxx[5])

    Vcur[5] += (dt) * (1/C) * (-gsyn45*sxx[6]*(Vprev[5]-Vsyn))
    sxx[6] += ds(Vprev[4], sxx[6])

    # Update the time
    t .= t .+ dt
end

function animstep!(Vcur, mNa, hNa, mCaT, hCaT, mCaS, hCaS, mA, hA,
                   mKCa, mKd, mH, Ca, gCaS, gCaSi, zCaS, gA, gAi, zA, gCaS_r, gA_r,
                   Iapp, gNa, gCaT, gKCa, gKd, gH, gleak, gsyn21, gsyn12,
                   gsyn13, gsyn53, gsyn54, gsyn45, gEl23, gEl43, tmKCa,
                   sxx, C, α, β, Kp, Ki, Kt, gsth, gsth_prev, guth, Vth, u_max, t,
                   trace1, trace2, trace3, trace4, trace5, tail)
    # Compute 100 points to be plotted to notify trace (accelerate ploting time)
    for j = 1 : 100
        # Only plot 1 computed point in 0.2ms
        for i = 1 : Int64(0.2/dt)
            progress_for_one_step!(Vcur, mNa, hNa, mCaT, hCaT, mCaS, hCaS, mA, hA,
                                   mKCa, mKd, mH, Ca, gCaS, gCaSi, zCaS, gA, gAi, zA, gCaS_r, gA_r,
                                   Iapp, gNa, gCaT, gKCa, gKd, gH, gleak, gsyn21, gsyn12,
                                   gsyn13, gsyn53, gsyn54, gsyn45, gEl23, gEl43, tmKCa,
                                   sxx, C, α, β, Kp, Ki, Kt, gsth, gsth_prev, guth, Vth, u_max, t)
        end
        # If we reached the end of the box, just erase the figure
        if abs((t[1]%5000)/1e3 - 5) < 0.03
            for k = 1 : tail + 1
                push!(trace1[], Point2f(0., -100.))
                push!(trace2[], Point2f(0., -100.))
                push!(trace3[], Point2f(0., -100.))
                push!(trace4[], Point2f(0., -100.))
                push!(trace5[], Point2f(0., -100.))
            end
        else # add the computed point to the trace
            push!(trace1[], Point2f((t[1]%5000)/1e3, Vcur[1]))
            push!(trace2[], Point2f((t[1]%5000)/1e3, Vcur[2]))
            push!(trace3[], Point2f((t[1]%5000)/1e3, Vcur[3]))
            push!(trace4[], Point2f((t[1]%5000)/1e3, Vcur[4]))
            push!(trace5[], Point2f((t[1]%5000)/1e3, Vcur[5]))
        end
    end
    # Notify to update the plot
    notify(trace1)
    notify(trace2)
    notify(trace3)
    notify(trace4)
    notify(trace5)
end

# Initial conditions
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
mKCa = mKCa_inf.(Vprev, Ca)
mKd = mKd_inf.(Vprev)
mH = mH_inf.(Vprev)

gCaS_r = zeros(5)
gA_r = zeros(5)
for j = 1 : 5
    (gCaS_r[j], gA_r[j]) = DICs_gmax_neuromodCaSA(gNa[j], gCaT[j], gKd[j], gKCa[j], gH[j],
                                                  gleak[j], gs_th[j], gu_th[j], Vth[j], tmKCa[j])
end

s21 = s_inf(Vprev[2])
s12 = s_inf(Vprev[1])
s13 = s_inf(Vprev[1])
s53 = s_inf(Vprev[5])
s54 = s_inf(Vprev[5])
s45 = s_inf(Vprev[4])
sxx = [s21, s12, s13, s53, s54, s45]

# Initialize time
t = zeros(1)
t .= 0.

# Initialize the Makie figure
(tail, trace1, trace2, trace3, trace4, trace5, fig, screen) = makefig(Vcur, t)

# Initialize the cursor
lsgrid = SliderGrid(fig[6, 5:30], (label="", range=-8:0.1:5,
                    startvalue=-8., format=""))

# Add listener
gsth = lsgrid.sliders[1].value
gs_th_prev = [gsth_sim, gsth_sim, gsth_sim, gsth[], gsth[]]

# Creating labels for the slider so that fontsize is tunable (not supported by SliderGrid)
txtTest = Label(fig[6, 3], text=L"g_s(V_{th}) = ", tellwidth=false, fontsize=50, halign=:right)
str_gsth = lift(gs -> "$gs", gsth)
txtTest = Label(fig[6, 4], text=str_gsth, tellwidth=false, fontsize=35, halign=:left)

# Creating the button
run = Button(fig[7, 1:30]; label="Run/stop", tellwidth=false, height=60, width=250, fontsize=30)

# Creating an observable flag to know whether simulation is running or not
isrunning = Observable(false)

# Add a first button instruction to toggle flag value
on(run.clicks) do clicks; isrunning[] = !isrunning[]; end

# Add a second button instruction to launch the simulation when flag is true and figure opened
on(run.clicks) do clicks
    @async while isrunning[]
        isopen(fig.scene) || break
        animstep!(Vcur, mNa, hNa, mCaT, hCaT, mCaS, hCaS, mA, hA,
                  mKCa, mKd, mH, Ca, gCaS, gCaSi, zCaS, gA, gAi, zA, gCaS_r, gA_r,
                  Iapp, gNa, gCaT, gKCa, gKd, gH, gleak, gsyn21, gsyn12,
                  gsyn13, gsyn53, gsyn54, gsyn45, gEl23, gEl43, tmKCa,
                  sxx, C, α, β, Kp, Ki, Kt, [gsth_sim, gsth_sim, gsth_sim, gsth[], gsth[]],
                  gs_th_prev, gu_th, Vth, u_max, t, trace1, trace2, trace3, trace4, trace5, tail)
        sleep(0.00000001)
    end
end
