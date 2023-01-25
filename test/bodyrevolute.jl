using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra

include(joinpath(@__DIR__, "components.jl"))

@named revol = Revolute()
@named bod = Body()
@named bod2 = Body()
@named bod3 = Body()
@named fix = FixedTranslation()
@named fix2 = FixedTranslation()
@named pos = Position()

ges_eqs = [
      # connect(bod.fa, fix2.fa)
      # connect(revol.fa,fix2.fb)
      # connect(revol.fb, fix.fa)
      # connect(fix.fb,bod2.fa)
      # connect(revol.fla,pos.flb)
      connect(bod.fa, fix2.fa)
      connect(fix.fb,fix.fa)
      connect(fix.fb, bod3.fa)
      connect(revol.fa, fix.fb)
      connect(revol.fb,bod2.fa)
      connect(revol.fla,pos.flb)
  ]

  #bod fix2 revol fix bod reihenfolge

@named _vehicle =ODESystem(ges_eqs, t)
@named vehicle =compose(_vehicle,[revol,bod,bod3,bod2,pos,fix,fix2])
sys = structural_simplify(vehicle)



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


