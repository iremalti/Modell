using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra

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
@named angle = AngleSensor()
@named pos = Position()
@named trap = Trapezoid()

ges_eqs = [
            connect(inertia.flb, rad.fla)
            connect(inertia2.flb, rad2.fla)
            connect(tor.flb, inertia2.fla)
            connect(bod.fa, fix2.fb)
            connect(fix2.fb, fix.fa)
            connect(fix.fb, revol.fa)
            connect(rad2.fa, fix2.fa)
            connect(revol.fb, rad.fa)
            connect(lastverteilung.ff, rad.dl) 
            connect(lastverteilung.fr, rad2.dl)
            connect(lastverteilung.frl, rad2.flo) 
            connect(lastverteilung.fl, rad.flo) 
            connect(lastverteilung.flat, rad.fat)  
            connect(konst.y, tor.in)  
            connect(konst.y, lastverteilung.M)  
            connect(angle.out,lastverteilung.d)
            connect(angle.fla,pos.flb)
            connect(pos.flb,revol.fla)          
            connect(trap.x, pos.phir)
          ]
# flang wird für revol  verwendet um t= 0 zu vermeiden
           

@named _vehicle =ODESystem(ges_eqs, t)
@named vehicle =compose(_vehicle,[fix,fix2,tor,lastverteilung,rad,rad2,revol,bod,inertia,inertia2,konst,angle,pos,trap,])
sys = structural_simplify(vehicle)

include(joinpath(@__DIR__, "printbalance.jl"))

u0 = [    ]
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Tsit5())
plot(sol)








###Irrelevante Gleichungen

#   f =fN * u_F;                                                                # irrelevant?? wird nicht weiter benutzt
#   lossPower = f*v_slip;                                                       # irrelevant?? wird nicht weiter benutzt
#   vAdhesion = noEvent(max(sAdhesion*abs(radius*w_roll),vAdhesion_min));       # irrelevant?? wird nicht weiter benutzt

#   vSlide = noEvent(max(sSlide*abs(radius*w_roll),vSlide_min));                # irrelevant?? wird nicht weiter benutzt

