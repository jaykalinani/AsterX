module BenchmarkPlots

using AbstractPlotting
using CSV
using DataFrames
using Makie



# Load data

iterations = 10
integrator_steps = 4

bm = DataFrame(CSV.File("benchmark-minkowski-unigrid.csv"))

bm.points = bm.size.^3
bm.cores = bm.nodes .* bm.cores_node
bm.cost = bm.time .* bm.cores
bm.pointcost = bm.cost ./ bm.points ./ (iterations * integrator_steps)



# Plot data

sc = Scene()
scatter!(sc, log10.(bm.cores), log10.(bm.pointcost))
xlims!(sc,
       (floor(minimum(log10.(bm.cores))),
        ceil(maximum(log10.(bm.cores)))))
xlabel!(sc, "log₁₀(cores)")
ylims!(sc,
       (floor(minimum(log10.(bm.pointcost))),
        ceil(maximum(log10.(bm.pointcost)))))
ylabel!(sc, "log₁₀(seconds per grid point)")

end
