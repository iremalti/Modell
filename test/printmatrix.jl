
ts = TearingState(expand_connections(vehicle))
im = ModelingToolkit.incidence_matrix(ts.structure.graph, 1.0f0)
using ModelingToolkit.BipartiteGraphs
ts = TearingState(expand_connections(rad))
inc_org = BipartiteGraphs.incidence_matrix(ts.structure.graph)
using GLMakie
function show_incidence_matrix(IM)
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        aspect = DataAspect(),
        ylabel="equation id",
        xlabel="state/unknown/variable id",
        xgridcolor=:white,
        ygridcolor=:white,
        xticks=1:1:size(IM, 2),
        yticks=(1:1:size(IM, 1), string.(-collect(1:1:size(IM, 1)).+size(IM, 1))) 
        )

    ima = image!(ax, Matrix(rotr90(IM)), interpolate=false)
    GLMakie.translate!(ima, 0, 0, -100)

    display(fig)
end
show_incidence_matrix(inc_org)




foreach(equations(expand_connections(vehicle))) do eq; println(eq) end