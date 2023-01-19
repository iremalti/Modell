using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
#using Symbolics

@variables t
D= Differential(t)

@register_symbolic smooth_pole(w_roll, grenze)
function smooth_pole(w_roll, grenze )
    if (abs(w_roll) < grenze) 
        y = (-1/(2 * abs(grenze)^3) * w_roll^2 + 3/(2 * abs(grenze)));
        else
        y= (1/abs(w_roll));
        end
    return y
end

#Connectoren entsprechend der Modelica Komponenten

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


function FixedTranslation(;name, r =[0. , 0.5], r0 =[0. , 0.],_R=[0. 0.;0. 0.], k =[0. , 0.])     
    @named fa= frame_a() 
    @named fb= frame_b()
    ps= @parameters rx(t)=r[1] ry(t)=r[2]
    sts= @variables  r0x(t)=r0[1] r0y(t)=r0[2] Rx(t)= _R[1] Rp(t)= _R[2] Ry(t)= _R[3] Rq(t)= _R[4] kx(t)=k[1]  ky(t)=k[2]      
    eqs= [  Rx ~ cos(fa.phi)
            Rp ~ -sin(fa.phi)
            Ry ~ sin(fa.phi)
            Rq ~ cos(fa.phi)
            kx ~ fb.fy
            ky ~ -fb.fx
            Rx * rx + Rp * ry .~ r0x                                                       # Position
            Ry * rx + Rq * ry .~ r0y                                                
            fa.ox .+ r0 .~ fb.ox                                                           # Positionen verbinden
            fa.oy .+ r0 .~ fb.oy
            fa.phi .~ fb.phi                                                               # winkel verbinden
            fa.fx .+ fb.fx .~ 0                                                            # kräfte ausgleich
            fa.fy .+ fb.fy .~ 0
            fa.tau .+ fb.tau .+ r0x * kx .+ r0y * ky .~ 0                                  # Moment ausgleich
            ]
    compose(ODESystem(eqs, t, sts, ps; name=name), [fa,fb])
end


function Revolute(;name, phi=0.0, w=0.0, z=0.0, tau=0.0, o=[0.,0.])     
    @named fa= frame_a() 
    @named fb= frame_b()
    @named fla= flange_a()
    sts= @variables phi(t) w(t) z(t) tau(t)  ox(t)=o[1] oy(t)=o[2]
    eqs= [  D.(phi) ~ w                                                                    # Differentialgleichungen
            D.(w) ~ z 
            ox ~ fa.ox
            oy ~ fa.oy
            ox ~ fb.ox
            oy ~ fb.oy
            fa.ox ~ fb.ox                                                                  # anfangs und endconnector bindungen
            fa.oy ~ fb.oy 
            fa.phi + phi ~ fb.phi                                                          # phi erneuert an den connector übergeben
            fa.fx + fb.fx ~ 0                                                              # kräfteausgleich
            fa.fy + fb.fy ~ 0 
            fa.tau + fb.tau ~ 0                                                            # moment ausgleich
            fa.tau ~ tau                                                                   # moment festlegen
            ]                
    compose(ODESystem(eqs, t, sts, []; name=name), [fa,fb,fla])
end
 

function Rad(;name, _R= [0. 0. ;0. 0.], v=[0., 0.], phi_roll=0.0, w_roll=0.0, e0=[0.,0.], v_lat=0.0, v_long=0.0, v_slip_lat=0.0, v_slip_long=0.0,v_slip=0.0,fN=0.0, f_long=0.0, f_lat=0.0, s_long=0.0, s_lat=0.0, s_F=0.0,u_F=0.0, u_long=0.0,u_lat=0.0)
    @named fa= frame_a() 
    @named fla= flange_a()
    @named dl = realInput()
    @named flo = realOutput()
    @named fat= realOutput()

    ps= @parameters radius=0.3 C=1.6 B=7 K=1 grenze=0.01   ex(t)=1 ey(t)=1  # e=[1,1]
    sts= @variables Rx(t)= _R[1] Rp(t)= _R[2] Ry(t)= _R[3] Rq(t)= _R[4] vx(t)=v[1] vy(t)=v[2] phi_roll(t) w_roll(t)  e0x(t)=e0[1] e0y(t)=e0[2] v_lat(t) v_long(t) v_slip_lat(t) v_slip_long(t) v_slip(t) fN(t) f_long(t) f_lat(t) s_long(t) s_lat(t) s_F(t) u_F(t) u_long(t) u_lat(t)
    eqs= [  Rx ~ cos(fa.phi)
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


function Lastenverteilung(;name)     
   
    @named fl= realInput()
    @named flat= realInput()
    @named d= realInput()
    @named frl= realInput()
    @named ff= realOutput()
    @named fr= realOutput()
    @named M= realInput()
 
    ps= @parameters  m = 1450 g=9.81 h=0.4 lr=1.59 lf=1.1  
    eqs= [ M.x ~ ff.x + fr.x - m * g                                                                                    # Gleichgewicht der Momente
           fr.x ~ ( h *(fl.x * cos(d.x) - flat.x * sin(d.x)+ frl.x) + ff.x * lf ) / lr         # Gleichgewicht der Kräfte
            ]                
    compose(ODESystem(eqs, t, [], ps; name=name), [fl,flat,d,frl,ff,fr,M])
end

# entsprechen den Modelica Komponenten

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

function Torque(;name)   
    @named in =realInput()
    @named flb = flange_b()
    sts = @variables tau(t)
    eqs = [tau ~ in.x
           -in.x ~ flb.tau ]
    compose(ODESystem(eqs, t, sts,[]; name=name), [flb,in])
end

function Konst(;name)   
    @named y = realOutput()
    ps = @parameters k = 100
    eqs = [
            y.x ~ k
          ]
    compose(ODESystem(eqs, t, [],ps; name=name), y)
end

# function ende(;name)   
#     @named y = realInput()
#     ps = @parameters k = 0
#     eqs = [
#             y.x ~ k
#           ]
#     compose(ODESystem(eqs, t, [],ps; name=name), y)
# end

function AngleSensor(;name)   
    @named out = realOutput()
    @named fla = flange_a()
    eqs = [ 
           out.x ~ fla.phi 
           fla.tau ~ 0
          ]
    compose(ODESystem(eqs, t, [],[]; name=name), [fla,out])
end

function Position(;name, phi=0.0, w=0.0, z=0.0, w_crit=0.0)   
    @named phir =realInput()
    @named flb = flange_b()
    ps= @parameters f_crit=50  af=1.3617  bf=0.6180
    sts = @variables phi(t) w(t) z(t) w_crit(t) 
    eqs = [ w_crit ~ pi * f_crit
            phi ~ phir.x
            phi ~ flb.phi
            D.(phi) ~ w                                                                    # Differentialgleichungen
            D.(w) ~ z 
            z ~ ((phir.x - phi)* w_crit - af * w) * (w_crit/bf)
          ]
    compose(ODESystem(eqs, t, sts, ps; name=name), [flb,phir])
end

function Trapezoid(;name,offset =pi/4)#,amplitude, rising, width, falling, period, nperiod, startTime)   
    @named x =realOutput()
   
    ps= @parameters offset=offset amplitude=0*pi/4 rising=0.1 width=0.8 falling=0.1 period=1.8 nperiod=-1 startTime=0
    eqs = [ 
             x.x ~ offset + signal(amplitude, rising, width, falling, period, nperiod, startTime)
          ]
    compose(ODESystem(eqs, t, [],ps; name=name), x)
end

#offset braucht ein default



@register_symbolic zeit(t)
value_vector = LinRange(0., 10., 10)     
zeit(t) = t >=10 ? alue_vector[end] : value_vector[Int(floor(t))+1]


@register_symbolic signal( amplitude, rising, width, falling, period, nperiod, startTime)
function signal(amplitude, rising, width, falling, period, nperiod, startTime)
    T_width = rising + width
    T_falling = T_width + falling
    count = integer((time - startTime)/period)
    T_start = startTime + count * period
    if (zeit(t) < startTime || nperiod == 0 || (nperiod > 0 && count >= nperiod)) return 0 
    elseif (zeit(t)< T_start + T_rising) return amplitude*(time - T_start)/rising 
    elseif (zeit(t) < T_start + T_width) return amplitude 
    elseif (zeit(t) < T_start + T_falling) return amplitude*(T_start + T_falling - time)/falling 
    else return 0
    end
end

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
#@named e = ende()

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
            connect(pos.flb,revol.fla)          # flang wird für revol  verwendet um t= 0 zu vermeiden
            connect(trap.x, pos.phir)
          ]


@named _vehicle =ODESystem(ges_eqs, t)
@named vehicle =compose(_vehicle,[fix,fix2,tor,lastverteilung,rad,rad2,revol,bod,inertia,inertia2,konst,angle,pos,trap,])



function print_balance(model)
    println(string(model.name))
    println("\tequations: $(length(equations(expand_connections(model))))")
    println("\tstates: $(length(states(model)))")
    println("-"^64)
    for s in model.systems
        neq = length(equations(expand_connections(s)))
        nst = length(states(s))
        println("\t$(string(s.name))")
        println("\t\tequations: $neq")
        println("\t\tstates: $nst")
    end
end

print_balance(vehicle)

sys = structural_simplify(vehicle)

u0 = [    ]
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Tsit5())
plot(sol)


###Irrelevante Gleichungen

#   f =fN * u_F;                                                                # irrelevant?? wird nicht weiter benutzt
#   lossPower = f*v_slip;                                                       # irrelevant?? wird nicht weiter benutzt
#   vAdhesion = noEvent(max(sAdhesion*abs(radius*w_roll),vAdhesion_min));       # irrelevant?? wird nicht weiter benutzt

#   vSlide = noEvent(max(sSlide*abs(radius*w_roll),vSlide_min));                # irrelevant?? wird nicht weiter benutzt

