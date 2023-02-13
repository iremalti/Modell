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