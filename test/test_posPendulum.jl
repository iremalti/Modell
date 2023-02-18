using ModelingToolkit, DifferentialEquations, LinearAlgebra
using GLMakie

include(joinpath(@__DIR__, "components.jl"))
include(joinpath(@__DIR__, "utils.jl"))

##

@named fixed = Fixed()
@named revolute = Revolute()
@named trans = FixedTranslation()
@named body = Body()
@named pos = Position()
#@named torque = Torque()
@named konst = Konst(k=pi)

eqs = [
      connect(fixed.fa, revolute.fa),
      connect(revolute.fb, trans.fa),
      connect(trans.fb, body.fa),
      connect(revolute.fla, pos.flb),
      connect(pos.phir, konst.y)
      #connect(torque.flb, revolute.fla),
      #connect(konst.y, torque.in)
]


@named _model =ODESystem(eqs, t)
@named model =compose(
    _model,
    [
        fixed,
        revolute,
        trans,
        body,
        pos,
        #torque,
        konst
    ]
)

modelPos = model

sys = structural_simplify(model)


u0 = [
      revolute.phi => pi/2,
      D(revolute.phi) => 0.0
     ]
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())


fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, sol.t, [u[1] for u in sol.u])
display(fig)