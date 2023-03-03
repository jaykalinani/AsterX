module BenchmarkPlots

using CSV
using CairoMakie
using DataFrames
using Primes
using Printf
using SixelTerm

################################################################################

function regex_join(regexs, sep)
    res = r""
    for (n, regex) in enumerate(regexs)
        n > 1 && (res *= sep)
        res *= regex
    end
    return res
end

################################################################################

function get_machine()
    machine = readlines(`./simfactory/bin/sim whoami`)
    machine = machine[begin]
    return machine
end

function get_basedir(machine::AbstractString)
    mdb = readlines(`./simfactory/bin/sim print-mdb-entry $machine`)

    basedir = filter(line -> match(r"^ *basedir *=.*", line) ≢ nothing, mdb)
    basedir = basedir[begin]
    basedir = match(r"^ *basedir *= *(.*)", basedir).captures[1]

    allocation = filter(line -> match(r"^ *allocation *=.*", line) ≢ nothing, mdb)
    allocation = allocation[begin]
    allocation = match(r"^ *allocation *= *(.*)", allocation).captures[1]

    user = filter(line -> match(r"^ *user *=.*", line) ≢ nothing, mdb)
    user = user[begin]
    user = match(r"^ *user *= *(.*)", user).captures[1]

    basedir = replace(basedir, r"@ALLOCATION@" => allocation)
    basedir = replace(basedir, r"@USER@" => user)

    return basedir
end

################################################################################

version() = "2022-04-21T14:28"

function submit(machine::AbstractString)
    # benchmark = "z4c"
    benchmark = "grhydrox"

    iter = "0000"

    all_nlevels = [1]           # [1, 8]
    all_size = [128, 192]       # [128, 192, 256]
    all_nodes = [128]           # [1, 2, 4, 8, 16, 32]
    if machine ∈ ["symmetry", "symmetry-llvm"]
        all_threads = [40]          # per process
    elseif machine ∈ ["summit-gpu"]
        all_threads = [24]          # per process
    else
        @assert false
    end

    all_blocking_factor = [8]
    all_max_grid_size = [128]   # [32, 64, 128]
    all_max_tile_size_x = [1024000]
    if endswith(machine, "-gpu")
        all_max_tile_size_y = [16]
        all_max_tile_size_z = [32]
    else
        all_max_tile_size_y = [4]  # [4, 8, 16, 32, 64]
        all_max_tile_size_z = [2]  # [1, 2, 4, 8, 16, 32]
    end

    if endswith(machine, "-gpu")
        configuration = "sim-gpu"
    elseif endswith(machine, "-llvm")
        configuration = "sim-llvm"
    else
        configuration = "sim"
    end

    if machine in ["symmetry", "symmetry-llvm"]
        queue = "defq"          # "debugq", "defq", "preq"
    elseif machine in ["symmetry-gpu"]
        queue = "gpuq"
    elseif machine in ["summit", "summit-gpu"]
        queue = "batch"
    else
        @assert false
    end

    cmds = []

    for nlevels in all_nlevels, size in all_size
        for nodes in all_nodes, threads in all_threads
            size_x, size_y, size_z = let
                sizes = Int[size, size, size]
                for (prim, count) in reverse(factor(nodes).pe)
                    for iter in 1:count
                        val, dir = findmin(sizes)
                        sizes[dir] *= prim
                    end
                end
                @assert prod(sizes) == nodes * size^3
                sizes
            end

            if machine ∈ ["symmetry"]
                cores = 40                  # per node
                smt = 1                     # per core
            elseif machine ∈ ["summit-gpu"]
                cores = 36                  # per node
                smt = 4                     # per core
            else
                @assert false
            end
            procs = nodes * cores * smt # total

            for blocking_factor in all_blocking_factor,
                max_grid_size in all_max_grid_size,
                max_tile_size_x in all_max_tile_size_x,
                max_tile_size_y in all_max_tile_size_y,
                max_tile_size_z in all_max_tile_size_z

                name = join(["benchmark", benchmark, "i$iter", "l$nlevels", "sx$size_x", "sy$size_y", "sz$size_z", "n$nodes",
                             "c$cores", "p$procs", "t$threads", "bf$blocking_factor", "gs$max_grid_size", "tx$max_tile_size_x",
                             "ty$max_tile_size_y", "tz$max_tile_size_z"], "-")

                @info "Submitting $name..."
                run(pipeline(Cmd(`./simfactory/bin/sim purge $name`; ignorestatus=true); stderr=devnull, stdout=devnull))
                cmd = run(Cmd(["./simfactory/bin/sim", "--machine=$machine", "submit", "$name", "--configuration=$configuration",
                               "--parfile=repos/cactusamrex/BenchmarkPlots/par/benchmark-$benchmark.par",
                               "--replace=\$nlevels=$nlevels", "--replace=\$ncells_x=$size_x", "--replace=\$ncells_y=$size_y",
                               "--replace=\$ncells_z=$size_z", "--replace=\$blocking_factor=$blocking_factor",
                               "--replace=\$max_grid_size=$max_grid_size", "--replace=\$max_tile_size_x=$max_tile_size_x",
                               "--replace=\$max_tile_size_y=$max_tile_size_y", "--replace=\$max_tile_size_z=$max_tile_size_z",
                               "--ppn-used=$(cores*smt)", "--procs=$procs", "--num-threads=$threads", "--num-smt=$smt",
                               "--queue=$queue", "--walltime=1:0:0"]); wait=false)
                push!(cmds, cmd)
            end
        end
    end

    @info "Waiting for all submissions to finish..."
    wait.(cmds)

    @info "Done."
    return nothing
end

################################################################################

function analyse(machine::AbstractString)
    df = DataFrame(; machine=String[], time=Float64[], iter=String[], nlevels=Int[], size_x=Int[], size_y=Int[], size_z=Int[],
                   nodes=Int[], cores_node=Int[], threads=Int[], threads_process=Int[], blocking_factor=Int[], max_grid_size=Int[],
                   max_tile_size_x=Int[], max_tile_size_y=Int[], max_tile_size_z=Int[])
    # benchmark = "z4c"
    benchmark = "grhydrox"

    basedir = get_basedir(machine)
    for dir in readdir(basedir)
        m = match(regex_join([r"benchmark", benchmark, r"i(\d+)", r"l(\d+)", r"sx(\d+)", r"sy(\d+)", r"sz(\d+)", r"n(\d+)",
                              r"c(\d+)", r"p(\d+)", r"t(\d+)", r"bf(\d+)", r"gs(\d+)", r"tx(\d+)", r"ty(\d+)", r"tz(\d+)"], r"-"),
                  dir)
        m ≡ nothing && continue
        iter = m.captures[1]
        (nlevels, size_x, size_y, size_z, nodes, cores_node, threads, threads_process, blocking_factor, max_grid_size, max_tile_size_x, max_tile_size_y, max_tile_size_z) = map(str -> parse(Int,
                                                                                                                                                                                             str),
                                                                                                                                                                                m.captures[2:end])

        stdout = try
            readlines("$basedir/$dir/output-0000/stdout.txt")
        catch
            continue            # Simulation is not ready or failed
        end

        time_evolve = filter(line -> match(r" Evolve$", line) ≢ nothing, stdout)
        isempty(time_evolve) && continue # Simulation is not ready or failed
        time_evolve = time_evolve[end]
        time_evolve = split(time_evolve)
        @assert length(time_evolve) == 5
        time_evolve = time_evolve[2]
        time_evolve = parse(Float64, time_evolve)

        time_outputgh = filter(line -> match(r" OutputGH$", line) ≢ nothing, stdout)
        if isempty(time_outputgh)
            # Time is missing (probably too small)
            time_outputgh = 0.0
        else
            time_outputgh = time_outputgh[end]
            time_outputgh = split(time_outputgh)
            @assert length(time_outputgh) == 5
            time_outputgh = time_outputgh[2]
            time_outputgh = parse(Float64, time_outputgh)
        end

        time = time_evolve - time_outputgh

        push!(df,
              [machine, time, iter, nlevels, size_x, size_y, size_z, nodes, cores_node, threads, threads_process, blocking_factor,
               max_grid_size, max_tile_size_x, max_tile_size_y, max_tile_size_z])
    end

    CSV.write("benchmark-$benchmark-$machine.csv", df)

    return nothing
end

################################################################################

function plot(machine::AbstractString)
    benchmark = "z4c"
    # benchmark = "grhydrox"

    # Constants

    flop_per_point = 5462
    bytes_per_point = 2 * 21 * 8 # read/write, variables, double
    iterations = 10
    integrator_steps = 4

    # Load data

    if machine == "summit-gpu"
        # 2 flop/fma, 80 SM, 32 f64 cores/SM, clock, 6 V100/node
        flop_per_second_per_node = 2 * 80 * 32 * 1530e+6 * 6
        bytes_per_second_per_node = 900.0e+9 * 6
    end

    if machine == "symmetry"
        flop_per_second_per_node = 2 * 2 * 8 * 2.40e+9 * 40
    end
    if machine == "symmetry-gpu"
        flop_per_second_per_gpu = 4 * 68 * 1545e+6 / 4
    end

    bm = DataFrame(CSV.File("data/benchmark-$benchmark-$machine.csv"))
    filter!(row -> row.time ≢ missing, bm)

    bm.work = flop_per_point * bm.nlevels .* bm.size_x .* bm.size_y .* bm.size_z * integrator_steps * iterations
    bm.cost = bm.nodes * flop_per_second_per_node .* bm.time
    # bm.work = bytes_per_point * bm.nlevels .* bm.size_x .* bm.size_y .* bm.size_z * integrator_steps * iterations
    # bm.cost = bm.nodes * bytes_per_second_per_node .* bm.time
    bm.efficiency = bm.work ./ bm.cost
    bm.load = bm.work ./ bm.nodes

    # Plot data

    fig = Figure(; resolution=(1024, 576))
    axis = Axis(fig[1, 1]; title="Scaling benchmark ($machine)", xlabel="nodes", xscale=log10,
                xtickformat=xs -> [(@sprintf "%.0f" x) for x in xs],
                ylabel="Efficiency (algorithmic flop / hardware flop)",
                # ylabel="Efficiency (algorithmic bytes / hardware bytes)",
                ytickformat=ys -> [(@sprintf "%.1f%%" (100 * y)) for y in ys])
    machine_color = Dict("symmetry" => :red, "symmetry-llvm" => :green, "symmetry-gpu" => :blue, "summit" => :pink,
                         "summit-gpu" => :orange)
    xmin = 1
    xmax = maximum(bm.nodes)
    ymin = 0.0
    ymax = maximum(bm.efficiency)

    # Plot all data points
    scatter!(axis, bm.nodes, bm.efficiency; color=machine_color[machine])

    # # Connect highest efficiency
    # let
    #     xvals = Int[]
    #     yvals = Float64[]
    #     for nodes in sort!(unique(bm.nodes))
    #         bm_node = filter(row -> row.nodes == nodes, bm)
    #         efficiency = maximum(bm_node.efficiency)
    #         push!(xvals, nodes)
    #         push!(yvals, efficiency)
    #     end
    #     scatterlines!(axis, xvals, yvals; color=machine_color[machine], linestyle=:dash)
    # end

    # Connect strong scaling points (those with equal work)
    for work in unique(bm.work)
        bm_work = filter(row -> row.work == work, bm)
        xvals = Int[]
        yvals = Float64[]
        for nodes in sort!(unique(bm_work.nodes))
            bm_work_node = filter(row -> row.nodes == nodes, bm_work)
            efficiency = maximum(bm_work_node.efficiency)
            push!(xvals, nodes)
            push!(yvals, efficiency)
        end
        scatterlines!(axis, xvals, yvals; color=machine_color[machine], linestyle=:dash)
    end

    # Connect weak scaling points (those with equal lead, i.e. work per node)
    for load in unique(bm.load)
        bm_load = filter(row -> row.load == load, bm)
        xvals = Int[]
        yvals = Float64[]
        for nodes in sort!(unique(bm_load.nodes))
            bm_load_node = filter(row -> row.nodes == nodes, bm_load)
            efficiency = maximum(bm_load_node.efficiency)
            push!(xvals, nodes)
            push!(yvals, efficiency)
        end
        scatterlines!(axis, xvals, yvals; color=machine_color[machine])
    end

    # Legend:
    # 0: GCC
    # 1: LLVM
    # 2: GPU

    # xlims!(1, 1.1 * xmax)
    ylims!(ymin, 1.1 * ymax)

    save("figures/benchmarks-$machine.png", fig)

    return fig
end

end
