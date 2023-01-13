using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics

@variables t
D= Differential(t)

@connector function frame_a(;name, pos=[0., 0.], phi=0.0,  _f=[0., 0.], tau=0.0)
    sts = @variables ox(t)=pos[1] oy(t)=pos[2] phi(t) fx(t)=_f [connect=Flow] fy(t)=_f [connect=Flow] tau(t) [connect=Flow]  
    ODESystem(Equation[], t, sts, []; name=name)
end

@connector function frame_b(;name, pos=[0., 0.], phi=0.0,  _f=[0., 0.], tau=0.0)
   sts = @variables ox(t)=pos[1] oy(t)=pos[2] phi(t) fx(t)=_f [connect=Flow] fy(t)=_f [connect=Flow] tau(t) [connect=Flow]  
   ODESystem(Equation[], t, sts, []; name=name)
end

@connector function flange_a(;name)
    sts = @variables phi(t)  tau(t) [connect=Flow]  
    ODESystem(Equation[], t, sts, []; name=name, defaults = Dict(phi => 0.0, tau => 0.0))
end

@connector function flange_b(;name)
    sts = @variables phi(t)  tau(t) [connect=Flow]  
    ODESystem(Equation[], t, sts, []; name=name, defaults = Dict(phi => 0.0, tau => 0.0))
end

@connector function realInput(;name, x=0.0)  
    sts = @variables x(t) 
    ODESystem(Equation[], t, sts, []; name=name)
end

@connector function realOutput(;name, x=0.0)  
    sts = @variables x(t)  
    ODESystem(Equation[], t, sts, []; name=name)
end

@register_symbolic smooth_pole(w_roll, grenze)
function smooth_pole(w_roll, grenze )
    if (abs(w_roll) < grenze) 
        y = (-1/(2 * abs(grenze)^3) * w_roll^2 + 3/(2 * abs(grenze)));
        else
        y= (1/abs(w_roll));
        end
    return y
end

function Body(;name, v=[0., 0.], a=[0., 0.],phi=0.0, w=0.0, z=0.0)  
    @named fa = frame_a()
    ps = @parameters m = 1450 I=1.8 
    sts= @variables vx(t)=v[1] vy(t)=v[2]  ax(t)=v[1] ay(t)=a[2]  phi(t) w(t) z(t)
    
    eqs = [fa.tau ~ I * z                                                                 # Drallsatz
           D.(fa.ox) .~ vx                                                                # Geschwindigkeitsableitung 
           D.(fa.oy) .~ vy                                                   
           D.(vx) .~ ax                                                                  # Beschleunigungsableitung   
           D.(vy) .~ ay                                                  
           D.(phi) .~ w                                                                  # Winkelgeschwindigkeit
           D.(w) .~ z                                                                     # winkelbeschleunigung
           fa.phi~ phi                                                                    # winkel anpassen
           fa.fx .~ m .* ax                                                               # newton's law
           fa.fy .~ m .* ay ]
    compose(ODESystem(eqs, t, sts, ps; name=name), fa)
    
end

function Inertia(;name,phi=0.0, w=0.0, z=0.0)
    @named fla = flange_a()
    @named flb = flange_b()
    ps =@parameters J =1.8
    sts = @variables phi(t) w(t) z(t)
    eqs = [
        phi ~ fla.phi;
        phi ~ flb.phi;
        D.(phi) ~ w                
        D.(w) ~ z 
        J * z ~ fla.tau + flb.tau;
          ]
    compose(ODESystem(eqs, t, sts, ps; name=name), [fla,flb])
end

function Rad(;name, _R= [0. 0. ;0. 0.], v=[0., 0.], phi_roll=0.0, w_roll=0.0, e0=[0.,0.], e=[1,1], v_lat=0.0, v_long=0.0, v_slip_lat=0.0, v_slip_long=0.0,v_slip=0.0,fN=0.0, f_long=0.0, f_lat=0.0, s_long=0.0, s_lat=0.0, s_F=0.0,u_F=0.0, u_long=0.0,u_lat=0.0)
    @named fa= frame_a() 
    @named fla= flange_a()
    @named dl = realInput()
    @named flo = realOutput()
    @named fat= realOutput()

    ps= @parameters radius=0.3 C=1.6 B=7 K=1 grenze=0.01
    sts= @variables Rx(t)= _R[1] Rp(t)= _R[2] Ry(t)= _R[3] Rq(t)= _R[4] vx(t)=v[1] vy(t)=v[2] phi_roll(t) w_roll(t) ex(t)=e[1] ey(t)=e[2] e0x(t)=e0[1] e0y(t)=e0[2] v_lat(t) v_long(t) v_slip_lat(t) v_slip_long(t) v_slip(t) fN(t) f_long(t) f_lat(t) s_long(t) s_lat(t) s_F(t) u_F(t) u_long(t) u_lat(t)

   # R= [cos(fa.phi), -sin(fa.phi) ,sin(fa.phi) ,cos(fa.phi)]

    eqs= [  dl.x ~ 100
            Rx ~ cos(fa.phi)
            Rp ~ -sin(fa.phi)
            Ry ~ sin(fa.phi)
            Rq ~ cos(fa.phi)
            D.(fa.ox) .~ vx                                                                      #Geschwindigkeitsableitung 
            D.(fa.oy) .~ vy 
            fla.phi ~ phi_roll
            D.(phi_roll) ~ w_roll   
            e0x ~ ex * Rx + ey * Rp                                                                 # einheitsvektor bestimmung
            e0y ~ ex * Ry + ey * Rq
            v_lat ~ -vx * e0y + vy * e0x
            v_long ~ vx * e0x + vy * e0y
            v_slip_lat ~ v_lat
            v_slip_long ~ v_long - radius * w_roll
            v_slip ~ sqrt(v_slip_long^2 + v_slip_lat^2)
            fN ~ dl.x                                                                   # Normalkraft Bestimmung
            f_long ~ fa.fx * e0x + fa.fy * e0y         
            f_lat ~ fa.fy * e0x - fa.fx * e0y
            -f_long * radius ~ fla.tau   
            fa.tau ~ 0         
            s_F ~ sqrt(s_long^2 + s_lat^2)
            u_F ~ K *sin(C* atan(B*s_F))                                                             # magic formular
            u_long ~ (-s_long/s_F) * u_F                                                             # Gesamtreibungskräfte
            u_lat ~ (-s_lat/s_F) * u_F           
            s_long ~ -(v_slip_long /radius) * smooth_pole(w_roll, grenze)    
            s_lat ~ -(v_slip_lat /radius) * smooth_pole(w_roll, grenze)
            f_long ~ u_long * fN                                                                     # kräfte neu berechnen
            f_lat ~ u_lat * fN
            flo.x ~ f_long                                                                   #kräfte verbinden
            fat.x ~ f_lat
            ]           
    compose(ODESystem(eqs, t, sts, ps; name=name), [fa,fla,dl,flo,fat])
end


@named rad = Rad()
@named bod = Body()
@named inertia =Inertia()


ges_eqs = [
            connect(inertia.flb, rad.fla)
            connect(bod.fa, rad.fa)
          ]

@named _vehicle =ODESystem(ges_eqs, t)
@named vehicle =compose(_vehicle,[rad,bod,inertia])



subsystems = [
    rad,
    bod,
    inertia,
]
 
for ssys in subsystems
    println(string(nameof(ssys)))
    println("\tneqs = $(length(equations(expand_connections(ssys))))")
    println("\tnstates = $(length(states(ssys)))")
end

sys = structural_simplify(vehicle)
