using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics

@variables t
D = Differential(t)

@connector function frame_a(;name, pos=[0., 0.], phi=0.0, _f=[0., 0.], tau=0.0)
   
    sts0 = @variables pos(t)[1:2,1]=pos phi(t) f(t)[1:2,1]=_f [connect=Flow] tau(t) [connect=Flow]  
    sts1 =[Symbolics.scalarize.(sts0)...;]
    sts2 =[sts1...;]
    ODESystem(Equation[], t, sts2, []; name=name)
end

@connector function frame_b(;name, pos=[0., 0.], phi=0.0, _f=[0., 0.], tau=0.0)
    
    sts0 = @variables pos(t)[1:2,1]=pos phi(t) f(t)[1:2,1]=_f [connect=Flow] tau(t) [connect=Flow]  
    sts1 =[Symbolics.scalarize.(sts0)...;]
    sts2 =[sts1...;]
    ODESystem(Equation[], t, sts2, []; name=name)
end

# @connector function frame_b(;name, pos=[0., 0.], phi=0.0, _f=[0., 0.], tau=0.0)
#     sts = @variables pos(t)[1:2,1]=pos phi(t) f(t)[1:2,1]=_f [connect=Flow] tau(t) [connect=Flow]  
#     ODESystem(Equation[], t, [sts...;], []; name=name)
# end


# @connector function frame_a(;name, pos=[0., 0.], phi=0.0, fx=0.0, fy=0.0, tau=0.0)
#     sts = @variables ox(t)=pos[1] oy(t)=pos[2] phi(t)=phi fx(t)=fx [connect=Flow] fy(t)=fy [connect=Flow] tau(t)=tau [connect=Flow]  
#     ODESystem(Equation[], t, sts, []; name=name)
# end

# @connector function frame_b(;name, o=[0., 0.], phi=0.0, fx=0.0, fy=0.0, tau=0.0)
#     sts = @variables ox(t)=o[1] oy(t)=o[2] phi(t)=phi fx(t)=fx [connect=Flow] fy(t)=fy [connect=Flow] tau(t)=tau [connect=Flow]  
#     ODESystem(Equation[], t, sts, []; name=name)
# end



function body(;name, m = 1450, I=1.8, v=[0., 0.], a=[0., 0.],phi=0.0, w=0.0, z=0.0)  
    @named fa = frame_a()
    ps = @parameters m I 
    sts= @variables v(t)[1:2,1]=v  a(t)[1:2,1]=a phi(t) w(t) z(t)
    
    eqs = [fa.tau ~ I * z              #Drallsatz
           D.(fa.pos) .~ v             #Geschwindigkeitsableitung 
           D.(v) .~ a                  #Beschleunigungsableitung   
           D.(phi) .~ w                #Winkelgeschwindigkeit
           D.(w) .~ z                  #winkelbeschleunigung
           fa.phi~ phi                 #winkel anpassen
           fa.f .~ m .* a ]             #newton's law
    seqs = [Symbolics.scalarize.(eqs)...;]
    sts1 = [sts...;]
    sts2 = [sts1...;]
    compose(ODESystem(seqs, t, sts2, ps; name=name), fa)
    
end


function FixedTranslation(;name, r =[0. , 0.5], r0 =[0. , 0.],_R=[0. 0.;0. 0.], u =[0. , 0.])     
    @named fa= frame_a() 
    @named fb= frame_a()
    ps= @parameters r[1:2,1]= r
    sts= @variables  r0(t)[1:2,1]=r0 R(t)[1:2, 1:2]= _R  u(t)[1:2,1] =u
    R= [cos(fa.phi) -sin(fa.phi)
        sin(fa.phi) cos(fa.phi)]
    u = reshape([fb.f[2], -fb.f[1]], (1,2))
    eqs= [ R * r .~ r0                                                        
            fa.pos .+ r0 .~ fb.pos
            fa.phi .~ fb.phi
            fa.f .+ fb.f .~ 0
            fa.tau .+ fb.tau .+ u * r0 .~ 0]
    seqs = [Symbolics.scalarize.(eqs)...;]
    compose(ODESystem(seqs, t, sts, ps; name=name), n)
end

@named l = FixedTranslation()


    #global seqs, sts, sts1, sts2
   # sts1= [Symbolics.scalarize.(sts)...;]
   # sts1= [sts...;]
   # sts2= [sts1...;]