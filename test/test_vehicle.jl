using ModelingToolkit, DifferentialEquations, LinearAlgebra
using GLMakie

include(joinpath(@__DIR__, "components.jl"))
include(joinpath(@__DIR__, "utils.jl"))

##

@named radHinten = Rad()
@named radVorne = Rad()
@named revolute = Revolute()
@named body = Body()
@named inertiaHinten = Inertia()
@named inertiaVorne = Inertia()
@named verteilung = Lastenverteilung(m=1450.0, g=9.81, h=0.4, lr=1.59, lf=1.1)
@named torque = Torque()
@named stabHinten = FixedTranslation(r0=[verteilung[:lr], 0.0])
@named stabVorne = FixedTranslation(r0=[verteilung[:lf], 0.0])
@named winkelSensor = AngleSensor()
@named winkelVorgabe = Position()
#@named sollWinkel = Konst(k=pi/8)
@named sollWinkel = Trapezoid(offset=0.0, amplitude=pi/8, rising=1.0, width=1.0, falling=1.0, period=10.0, nperiod=-1.0, startTime=5.0)
@named wunschMoment = Konst(k=2.0)




ges_eqs = [
    # start at rear wheel (i.e. its actuated inertia) and follow the chain until revolute
    connect(wunschMoment.y, torque.in)
    connect(torque.flb, inertiaHinten.fla)
    connect(inertiaHinten.flb, radHinten.fla)
    connect(radHinten.fa, stabHinten.fa)
    connect(stabHinten.fb, stabVorne.fa)
    connect(stabVorne.fb, revolute.fa)

    # start at front wheel and follow the chain until revolute
    connect(inertiaVorne.flb, radVorne.fla)
    connect(radVorne.fa, revolute.fb)

    # add the main body
    connect(stabHinten.fb, body.fa)

    # add steering control (and a sensor)
    connect(sollWinkel.y, winkelVorgabe.phir)
    connect(winkelVorgabe.flb, revolute.fla)
    connect(winkelSensor.fla, winkelVorgabe.flb)

    # inputs to the normal load distribution
    connect(radVorne.flo, verteilung.f_longFront)
    connect(radVorne.fat, verteilung.f_latFront)
    connect(radHinten.flo, verteilung.f_longRear)
    connect(wunschMoment.y, verteilung.tau_rear)
    connect(winkelSensor.out, verteilung.delta_steer)

    # outputs from the normal load distribution
    connect(verteilung.f_normalFront, radVorne.dl)
    connect(verteilung.f_normalRear, radHinten.dl)           
]


@named _vehicle = ODESystem(ges_eqs, t)
@named vehicle = compose(_vehicle, [
   radHinten,
   radVorne,
   revolute,
   body,
   inertiaHinten,
   inertiaVorne,
   verteilung,
   torque,
   stabHinten,
   stabVorne,
   winkelSensor,
   winkelVorgabe,
   sollWinkel,
   wunschMoment
   ]
)
sys = structural_simplify(vehicle)

u0 = [
    radHinten.phi_roll => 0.0,
    inertiaHinten.w => 0.1,
    radVorne.phi_roll => 0.0,
    inertiaVorne.w => 0.1,
    verteilung.f_normalRear.x => 0.0
]

# params = [
#     sollWinkel.offset => 0.0,
#     sollWinkel.amplitude => pi/8,
#     sollWinkel.rising => 0.2,
#     sollWinkel.width => 0.8,
#     sollWinkel.falling => 0.2,
#     sollWinkel.period => 10.0,
#     sollWinkel.nperiod => -1.0,
#     sollWinkel.startTime => 5.0,
# ]

# append!(u0, params)

saveat = LinRange(0, 50, 2^9) |> collect

tstops = get_event_times(saveat, sollWinkel)

append!(saveat, tstop)

sort!(saveat)


prob = ODEProblem(sys, u0, (saveat[1], saveat[end]), discrete_events = [tstops => [body.vx ~ body.vx]]);
sol = solve(prob, Rodas4(), saveat=saveat, dtmax=0.05)

##

# https://docs.sciml.ai/ModelingToolkit/stable/basics/FAQ/#Getting-the-index-for-a-symbol
indexof(sym, syms) = findfirst(isequal(sym), syms)


fig = Figure()

gleft = fig[1, 1] = GridLayout()
gright = fig[1, 2] = GridLayout()

# plot COM xy trajectory
axtraj = Axis(gright[1, 1], ylabel="com_y (m)", xlabel="com_x (m)", aspect = DataAspect())
lines!(axtraj, sol[body.fa.ox], sol[body.fa.oy])

# ... and mark the points where events occur
idx_ox = indexof(body.fa.ox, states(sys))
idx_oy = indexof(body.fa.oy, states(sys))
scatter!(axtraj, [(u[idx_ox], u[idx_oy]) for u in sol(tstops).u])

# plot steering angle over time
axsteer = Axis(gleft[1, 1], ylabel="steering (rad)", xlabel="t (s)")
lines!(axsteer, sol.t, sol[sollWinkel.y.x])

# ... and mark the event times
vlines!(axsteer, tstop, linestyle=:dashdot)

# plot speed over time
axspeed = Axis(gleft[2, 1], ylabel="speed (km/h)", xlabel="t (s)")
lines!(axspeed,
    sol.t,
    sqrt.(sol[radHinten.vx].^2 + sol[radHinten.vy].^2) * 3.6
)

# plot torque over time
axtorque = Axis(gleft[3, 1], ylabel="torque (Nm)", xlabel="t (s)")
lines!(axtorque,
    sol.t,
    sol[wunschMoment.y.x]
)



#lines!(ax, sol.t, sol[verteilung.f_normalRear.x])

display(fig)
