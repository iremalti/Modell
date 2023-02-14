import Base: getindex

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


function systems(osys::ODESystem)
    for sys in osys.systems
        display(sys)
        println("-"^32)
    end
end


function getindex(sys::ODESystem, sym::Symbol)
    defaults_map = ModelingToolkit.get_defaults(sys)
    defaults_list = defaults_map |> keys |> collect

    psym = findfirst(defaults_list) do child
        string(child) == string(sym)
    end

    ModelingToolkit.get_defaults(sys)[defaults_list[psym]]
end


function get_event_times(saveat, trap::ODESystem)
    tstop = Float64[]

    period = trap[:period]
    nperiod = trap[:nperiod]

    if nperiod == 0 || trap[:amplitude] == 0.0 || (trap[:rising] + trap[:width] + trap[:falling]) == 0.0
        return tstop
    end

    T_start = trap[:startTime]
    T_risen = T_start + trap[:rising]
    T_flat = T_risen + trap[:width]
    T_fallen = T_flat + trap[:falling]
    
    if nperiod < 0
        # case: infinitely many periods
        eventTime_ub = saveat[end]
    else
        eventTime_ub = min(saveat[end], T_start + nperiod*period)
    end

    for T in (T_start, T_risen, T_flat, T_fallen)
        append!(
            tstop,
            max(T, saveat[1]):period:eventTime_ub
        )
    end

    return sort!(tstop)
end