@variables t
D= Differential(t)

#Connectoren entsprechend der Modelica Komponenten

@connector function frame_a(;name, pos=[0., 0.], phi=0.0,  _f=[0., 0.], tau=0.0)
     sts = @variables ox(t)=pos[1] oy(t)=pos[2] phi(t) fx(t)=_f[1] [connect=Flow] fy(t)=_f[2] [connect=Flow] tau(t) [connect=Flow]  
     ODESystem(Equation[], t, sts, []; name=name)
 end


@connector function frame_b(;name, pos=[0., 0.], phi=0.0,  _f=[0., 0.], tau=0.0)
    sts = @variables ox(t)=pos[1] oy(t)=pos[2] phi(t) fx(t)=_f[1] [connect=Flow] fy(t)=_f[2] [connect=Flow] tau(t) [connect=Flow]  
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
    sts = @variables x(t) [input=true]
    ODESystem(Equation[], t, sts, []; name=name)
end


@connector function realOutput(;name, x=0.0)  
    sts = @variables x(t) [output=true]
    ODESystem(Equation[], t, sts, []; name=name)
end


function Fixed(;name, pos=[0., 0.], phi=0.0)  
    @named fa = frame_a()
    ps = @parameters posx=pos[1] posy=pos[2] phi=phi
    
    eqs = [
        fa.ox ~ posx
        fa.oy ~ posy
        fa.phi ~ phi
    ]
    compose(ODESystem(eqs, t, [], ps; name=name), fa)
end


function Body(;name, v=[0., 0.], phi=0.0, w=0.0, z=0.0)  
    @named fa = frame_a()
    ps = @parameters m=1.0 I=1.0 
    sts= @variables vx(t)=v[1] vy(t)=v[2] ax(t) ay(t) phi(t)=0.0 [state_priority = 100] w(t)=w [state_priority = 100] z(t)=z
    
    eqs = [fa.tau ~ I * z                                                                 # Drallsatz
           D.(fa.ox) ~ vx                                                                # Geschwindigkeitsableitung 
           D.(fa.oy) ~ vy                                                   
           D.(vx) ~ ax                                                                  # Beschleunigungsableitung   
           D.(vy) ~ ay                                                  
           D.(phi) ~ w                                                                  # Winkelgeschwindigkeit
           D.(w) ~ z                                                                     # winkelbeschleunigung
           fa.phi ~ phi                                                                    # winkel anpassen
           fa.fx ~ m * ax                                                               # newton's law
           fa.fy ~ m * ay ]
    compose(ODESystem(eqs, t, sts, ps; name=name), fa)
end


function FixedTranslation(;name, r0 =[1. , 0.])     
    @named fa= frame_a() 
    @named fb= frame_b()
    ps= @parameters rx=r0[1] ry=r0[2] # Versatz von fb gegenüber fa zerlegt in Koordinaten von fa
    #sts= @variables  r0x(t)=r0[1] r0y(t)=r0[2] Rx(t)= _R[1] Rp(t)= _R[2] Ry(t)= _R[3] Rq(t)= _R[4]
    sts= @variables  r0x(t) r0y(t) Rx(t) Rp(t) Ry(t) Rq(t) 
    eqs= [  
        Rx ~ cos(fa.phi)
        Rp ~ -sin(fa.phi)
        Ry ~ sin(fa.phi)
        Rq ~ cos(fa.phi)
        Rx * rx + Rp * ry ~ r0x                                                       # Position
        Ry * rx + Rq * ry ~ r0y                                                
        fa.ox + r0x ~ fb.ox                                                           # Positionen verbinden
        fa.oy + r0y ~ fb.oy
        fa.phi ~ fb.phi                                                               # winkel verbinden
        fa.fx + fb.fx ~ 0                                                            # kräfte ausgleich
        fa.fy + fb.fy ~ 0
        # Momentenbilanz um den Ursprung von fa (berechnet in Koordinaten von fa)
        fa.tau + fb.tau + r0x * fb.fy - r0y * fb.fx ~ 0
        #fa.tau + fb.tau + rx * fb.fy - ry * fb.fx ~ 0
        ]
    compose(ODESystem(eqs, t, sts, ps; name=name), [fa,fb])
end


function Revolute(;name, phi=0.0, w=0.0, z=0.0, tau=0.0)     
    @named fa= frame_a() 
    @named fb= frame_b()
    @named fla= flange_a()
    sts= @variables phi(t)=phi [state_priority = 100] w(t)=w [state_priority = 100] z(t)=z tau(t)=tau  
    eqs= [  D.(phi) ~ w                                                                    # Differentialgleichungen
            D.(w) ~ z 
            fla.phi ~ phi
            fla.tau ~ tau
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
 

function Rad(;name, v=[0., 0.], phi_roll=0.0, w_roll=0.0, e=[1., 0.])
    @named fa= frame_a() 
    @named fla= flange_a()
    @named dl = realInput()
    @named flo = realOutput()
    @named fat = realOutput()

    _e = normalize(e)

    ps= @parameters radius=0.3 C=1.6 B=7 K=1 grenze=0.01 ex(t)=_e[1] ey(t)=_e[2]
    sts= @variables Rx(t) Rp(t) Ry(t) Rq(t) vx(t)=v[1] vy(t)=v[2] phi_roll(t) w_roll(t) e0x(t) e0y(t) v_lat(t) v_long(t) v_slip_lat(t) v_slip_long(t) v_slip(t) fN(t) f_long(t) f_lat(t) s_long(t) s_lat(t) s_F(t) u_F(t) u_long(t) u_lat(t)
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


function Lastenverteilung(;name, m, g, h, lr, lf)     
    @named f_longFront = realInput()
    @named f_latFront = realInput()
    @named delta_steer = realInput()
    @named f_longRear = realInput()
    @named f_normalFront = realOutput()
    @named f_normalRear = realOutput()
    @named tau_rear = realInput()
 
    ps = @parameters  m=m g=g h=h lr=lr lf=lf  
    eqs = [
        tau_rear.x ~ f_normalFront.x + f_normalRear.x - m * g                                                                                    # Gleichgewicht der Momente
        f_normalRear.x ~ ( h *(f_longFront.x * cos(delta_steer.x) - f_latFront.x * sin(delta_steer.x) + f_longRear.x) + f_normalFront.x * lf ) / lr         # Gleichgewicht der Kräfte
    ]                
    compose(
        ODESystem(eqs, t, [], ps; name=name),
        [f_longFront, f_latFront, delta_steer, f_longRear, f_normalFront, f_normalRear, tau_rear]
    )
end

# entsprechen den Modelica Komponenten

function Inertia(;name, phi=0.0, w=0.0, z=0.0)
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
    @named in = realInput()
    @named flb = flange_b()

    eqs = [-in.x ~ flb.tau,]
    
    compose(ODESystem(eqs, t, [],[]; name=name), [flb,in])
end


function Konst(;name, k=1.0)   
    @named y = realOutput()
    ps = @parameters k=k
    eqs = [
            y.x ~ k
          ]
    compose(ODESystem(eqs, t, [],ps; name=name), y)
end


function AngleSensor(;name)   
    @named out = realOutput()
    @named fla = flange_a()
    eqs = [ 
           out.x ~ fla.phi 
           fla.tau ~ 0
          ]
    compose(ODESystem(eqs, t, [],[]; name=name), [fla,out])
end


function Position(;name, phi=0.0, w=0.0, z=0.0)   
    @named phir = realInput()
    @named flb = flange_b()
    ps= @parameters f_crit=50  af=1.3617  bf=0.6180
    sts = @variables phi(t) w(t)=w z(t) w_crit(t) 
    eqs = [ w_crit ~ pi * f_crit
            #phi ~ phir.x
            phi ~ flb.phi
            D.(phi) ~ w                                                                    # Differentialgleichungen
            D.(w) ~ z 
            z ~ ((phir.x - phi)* w_crit - af * w) * (w_crit/bf)
          ]
    compose(ODESystem(eqs, t, sts, ps; name=name), [flb,phir])
end


function Trapezoid(;name, offset=0.0, amplitude=1.0, rising=0.1, width=0.2, falling=0.1, period=1.0, nperiod=-1, startTime=0.0)   
    @named y = realOutput()
   
    ps= @parameters offset=offset amplitude=amplitude rising=rising width=width falling=falling period=period nperiod=-1 startTime=0
    eqs = [ 
             y.x ~ trapezoid(t, offset, amplitude, rising, width, falling, period, nperiod, startTime)
          ]
    
    # TODO include appropriate events to fire between signal sections.
    # but which one is appropriate? https://docs.sciml.ai/ModelingToolkit/stable/basics/Events/
    compose(ODESystem(eqs, t, [], ps; name=name), y)
end


@register_symbolic trapezoid(t, offset, amplitude, rising, width, falling, period, nperiod, startTime)

function trapezoid(t, offset, amplitude, rising, width, falling, period, nperiod, startTime)
    T_width = rising + width
    T_falling = T_width + falling
    count = fld(t - startTime, period)
    T_start = startTime + count * period
    
    if t < startTime || nperiod == 0 || (nperiod > 0 && count >= nperiod)
        t >= 35.3
        return offset
    elseif t < T_start + rising
        return offset + amplitude*(t - T_start)/rising 
    elseif t < T_start + T_width
        return offset + amplitude 
    elseif t < T_start + T_falling
        return offset + amplitude*(T_start + T_falling - t) / falling 
    else
        t >= 35.3
        return offset
    end
end


@register_symbolic smooth_pole(w_roll, grenze)
function smooth_pole(w_roll, grenze)
    if (abs(w_roll) < grenze) 
        y = (-1/(2 * abs(grenze)^3) * w_roll^2 + 3/(2 * abs(grenze)));
        else
        y= (1/abs(w_roll));
        end
    return y
end
