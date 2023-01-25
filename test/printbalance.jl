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
