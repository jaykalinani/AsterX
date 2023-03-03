using BenchX
using CSV
using DataFrames

# TODO:
# 1. create and store set of benchmarks, then run them independently, and analyze them independently
# 2. put benchmark description CSV into simulation directory

physics = PhysicsParams(;
    name="z4c",
    configuration="sim-llvm",
    parfile="repos/cactusamrex/BenchX/par/benchmark-z4c.par",
    extra_physics_params=Dict(
        # "\$nlevels" => 8,
        # "\$ncells_x" => 128,
        # "\$ncells_y" => 128,
        # "\$ncells_z" => 128,
        "\$nlevels" => 1,
        "\$ncells_x" => 256,
        "\$ncells_y" => 256,
        "\$ncells_z" => 256,
    ),
)

run = RunParams(;
    machine="symmetry-llvm",
    queue="debugq",
    walltime_seconds=3600,
    nodes=1,
    cores_per_node=40,
    pus_per_core=1,
    processes=2,
    threads_per_process=20,
    smts_per_thread=1,
    extra_run_params=Dict(
        # "\$max_grid_size" => 32,
        "\$max_grid_size" => 128,
        "\$max_tile_size_x" => 1073741824,
        "\$max_tile_size_y" => 16,
        "\$max_tile_size_z" => 16,
    ),
)

# benchmark = Benchmark(physics, run)
# 
# status = get_run_status(benchmark)
# if status ≡ rs_unknown
#     submit_run(benchmark)
# end
# if status ≢ rs_finished
#     wait_for_run(benchmark)
# end
# 
# timing = read_run_timing(benchmark)

benchmarks = Benchmark[]
for max_tile_size_y in [2, 4] # [2, 4, 8, 16, 32]
    for max_tile_size_z in [2] # [2, 4, 8, 16, 32]
        if max_tile_size_y ≥ max_tile_size_z
            for threads_per_process in [20, 40] # [1, 2, 4, 5, 8, 10, 20, 40]
                run′ = RunParams(;
                    machine=run.machine,
                    queue=run.queue,
                    walltime_seconds=run.walltime_seconds,
                    nodes=run.nodes,
                    cores_per_node=run.cores_per_node,
                    pus_per_core=run.pus_per_core,
                    # processes=run.processes,
                    # threads_per_process=run.threads_per_process,
                    smts_per_thread=run.smts_per_thread,
                    processes=run.cores_per_node ÷ threads_per_process,
                    threads_per_process=threads_per_process,
                    extra_run_params=copy(run.extra_run_params),
                )
                run′.extra_run_params["\$max_tile_size_y"] = max_tile_size_y
                run′.extra_run_params["\$max_tile_size_z"] = max_tile_size_z

                benchmark = Benchmark(physics, run′)
                push!(benchmarks, benchmark)
            end
        end
    end
end

benchmarks′ = Benchmark[]
for benchmark in benchmarks
    status = get_run_status(benchmark)
    status ≡ rs_unknown && push!(benchmarks′, benchmark)
end
submit_runs(benchmarks′)

wait_for_runs(benchmarks)

results = BenchmarkResult[]
for benchmark in benchmarks
    timing = read_run_timing(benchmark)
    result = BenchmarkResult(benchmark, timing)
    println(results)
    push!(results, result)
end

sort!(results; by=result -> result.timing.cell_updates_per_second, rev=true)

dataframe = DataFrame(;
    # PhysicsParams
    name=AbstractString[],
    configuration=AbstractString[],
    config_id=AbstractString[],
    build_id=AbstractString[],
    parfile=AbstractString[],
    [Symbol(key) => typeof(val)[] for (key, val) in physics.extra_physics_params]...,
    # RunParams
    machine=AbstractString[],
    queue=AbstractString[],
    walltime_seconds=Float64[],
    nodes=Int[],
    cores_per_node=Int[],
    pus_per_core=Int[],
    processes=Int[],
    threads_per_process=Int[],
    smts_per_thread=Int[],
    [Symbol(key) => typeof(val)[] for (key, val) in run.extra_run_params]...,
    # Timing
    submitted=Bool[],
    success=Bool[],
    evolution_seconds=Float64[],
    evolution_compute_seconds=Float64[],
    evolution_output_seconds=Float64[],
    evolution_cell_updates=Float64[],
    evolution_iterations=Int[],
    cells=Float64[],
    cell_updates_per_second=Float64[],
)
for result in results
    dict = merge(
        convert(Dict{Symbol,Any}, result.benchmark.physics),
        convert(Dict{Symbol,Any}, result.benchmark.run),
        convert(Dict{Symbol,Any}, result.timing),
    )
    push!(dataframe, dict)
end
CSV.write("benchmarks.csv", dataframe)
