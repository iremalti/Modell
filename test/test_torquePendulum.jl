using ModelingToolkit, DifferentialEquations, LinearAlgebra
using GLMakie

include(joinpath(@__DIR__, "components.jl"))
include(joinpath(@__DIR__, "utils.jl"))

##

@named fixed = Fixed()
@named revolute = Revolute()
@named trans = FixedTranslation()
@named body = Body()
#@named pos = Position()
@named torque = Torque()
@named konst = Konst(k=1)

eqs = [
      connect(fixed.fa, revolute.fa),
      connect(revolute.fb, trans.fa),
      connect(trans.fb, body.fa),
      #connect(revolute.fla, pos.flb),
      #connect(pos.phir, konst.y)
      connect(torque.flb, revolute.fla),
      connect(konst.y, torque.in)
]

@named _model =ODESystem(eqs, t)
@named model =compose(
    _model,
    [
        fixed,
        revolute,
        trans,
        body,
        #pos,
        torque,
        konst
    ]
)

modelTorque = model

sys = structural_simplify(model, allow_parameter=false)

##

u0 = [
      # instantly crashing parametrization for trans
      trans.rx => 1.0,
      trans.ry => 0.0,

      # initially working parametrization for trans (until t ~ sqrt(2*pi), i.e. until the pendulum has rotated by pi/2)
      #trans.rx => 1/sqrt(2),
      #trans.ry => 1/sqrt(2),

      # supposedly sufficient initial conditions
      body.phi => 0.0,
      D(body.phi) => 0.0,
      D(D(body.phi)) => 0.0,

      # initial conditions required by MTK
      D(D(body.fa.oy)) => 1.0,
      D(D(trans.r0x)) => 0.0,
      trans.r0y => 0.0,
      D(trans.r0y) => 0.0,
]

prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4(), saveat=0:0.1:10)

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, sol.t, sol[body.phi])
display(fig)