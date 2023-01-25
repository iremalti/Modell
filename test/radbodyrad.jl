using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra,GLMakie

include(joinpath(@__DIR__, "components.jl"))


@named rad = Rad()
@named rad2 = Rad()
@named revol = Revolute()
@named bod = Body()
@named inertia =Inertia()
@named inertia2 =Inertia()
@named lastverteilung = Lastenverteilung()
@named tor = Torque()
@named fix = FixedTranslation()
@named fix2 = FixedTranslation()
@named konst = Konst()
@named konst2 = Konst()
@named angle = AngleSensor()
@named pos = Position()
@named trap = Trapezoid()

ges_eqs = [
            connect(inertia.flb, rad.fla)
            connect(inertia2.flb, rad2.fla)
            connect(bod.fa, fix2.fb)
            connect(fix2.fb, fix.fa)
            connect(fix.fb, revol.fa)
            connect(rad2.fa, fix2.fa)
            connect(revol.fb, rad.fa)
            connect(rad.dl,konst.y)
            connect(rad2.dl,konst2.y)
         #   connect(angle.fla,pos.flb)
          #  connect(pos.flb,revol.fla)          
           # connect(trap.x, pos.phir)
                    
          ]

          #rad fix fix2 revol rad Reihenfolge

@named _vehicle =ODESystem(ges_eqs, t)
@named vehicle =compose(_vehicle,[rad,rad2,bod,inertia,inertia2,fix,fix2,konst,konst2,revol])
sys = structural_simplify(vehicle)

