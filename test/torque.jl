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


@named inertia =Inertia()
@named tor = Torque()
@named konst = Konst()

ges_eqs = [
    connect(konst.y, tor.in)
    connect(tor.flb, inertia.fla)
  ]

@named _vehicle =ODESystem(ges_eqs, t)
@named vehicle =compose(_vehicle,[konst, inertia, tor,])


sys = structural_simplify(vehicle)
