using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra


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

function Body(;name, v=[0., 0.], a=[0., 0.],phi=0.0, w=0.0, z=0.0)  
   @named fa = frame_a()
   ps = @parameters m = 1450 I=1.8 
   sts= @variables vx(t)=v[1] vy(t)=v[2] ax(t)=v[1] ay(t)=a[2] phi(t) w(t) z(t)
   
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



@named fix = FixedTranslation()
@named bod = Body()
@named bod2 = Body()

ges_eqs = [
    connect(bod.fa, fix.fa)
    connect(fix.fb, bod2.fa)
  ]

@named _vehicle =ODESystem(ges_eqs, t)
@named vehicle =compose(_vehicle,[fix,bod,bod2,])


sys = structural_simplify(vehicle)
