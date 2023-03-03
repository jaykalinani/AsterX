module BenchX

using SHA
using YAML

################################################################################

export PhysicsParams
struct PhysicsParams
    name::AbstractString
    configuration::AbstractString
    config_id::AbstractString
    build_id::AbstractString
    parfile::AbstractString
    extra_physics_params::Dict{AbstractString,Union{AbstractString,Float64,Int}}
    function PhysicsParams(;
        name::AbstractString,
        configuration::AbstractString,
        config_id::Union{Nothing,AbstractString}=nothing,
        build_id::Union{Nothing,AbstractString}=nothing,
        parfile::AbstractString,
        extra_physics_params::Dict,
    )
        if config_id ≡ nothing
            cactusdir = find_cactusdir()
            config_id = chomp(read(joinpath(cactusdir, "configs", configuration, "CONFIG-ID"), String))
        end
        if build_id ≡ nothing
            cactusdir = find_cactusdir()
            build_id = chomp(read(joinpath(cactusdir, "configs", configuration, "BUILD-ID"), String))
        end
        return new(name, configuration, config_id, build_id, parfile, extra_physics_params)
    end
end
function Base.show(io::IO, physics::PhysicsParams)
    println(io, "PhysicsParams:")
    println(io, "    name: ", physics.name)
    println(io, "    configuration: ", physics.configuration)
    println(io, "    config_id: ", physics.config_id)
    println(io, "    build_id: ", physics.build_id)
    println(io, "    parfile: ", physics.parfile)
    println(io, "    extra_physics_params:")
    for (key, val) in sort!(collect(physics.extra_physics_params))
        println(io, "        ", key, ": ", val)
    end
    return nothing
end
function Base.convert(::Type{DICT}, physics::PhysicsParams) where {DICT<:AbstractDict}
    dict = DICT(
        :name => physics.name,
        :configuration => physics.configuration,
        :config_id => physics.config_id,
        :build_id => physics.build_id,
        :parfile => physics.parfile,
    )
    for (key, val) in physics.extra_physics_params
        dict[Symbol(key)] = val
    end
    return dict
end
function Base.hash(physics::PhysicsParams, h::UInt)
    return hash(
        physics.name,
        hash(
            physics.configuration,
            hash(physics.config_id, hash(physics.build_id, hash(physics.parfile, hash(physics.extra_physics_params, h)))),
        ),
    )
end

export RunParams
struct RunParams
    machine::AbstractString
    queue::AbstractString
    walltime_seconds::Float64
    nodes::Int
    cores_per_node::Int
    pus_per_core::Int
    processes::Int
    threads_per_process::Int
    smts_per_thread::Int
    extra_run_params::Dict{AbstractString,Union{AbstractString,Float64,Int}}
    function RunParams(;
        machine::AbstractString,
        queue::AbstractString,
        walltime_seconds::Union{AbstractFloat,Integer},
        nodes::Int,
        cores_per_node::Int,
        pus_per_core::Int,
        processes::Int,
        threads_per_process::Int,
        smts_per_thread::Int,
        extra_run_params::Dict,
    )
        return new(
            machine,
            queue,
            walltime_seconds,
            nodes,
            cores_per_node,
            pus_per_core,
            processes,
            threads_per_process,
            smts_per_thread,
            extra_run_params,
        )
    end
end
function Base.convert(::Type{DICT}, run::RunParams) where {DICT<:AbstractDict}
    dict = DICT(
        :machine => run.machine,
        :queue => run.queue,
        :walltime_seconds => run.walltime_seconds,
        :nodes => run.nodes,
        :cores_per_node => run.cores_per_node,
        :pus_per_core => run.pus_per_core,
        :processes => run.processes,
        :threads_per_process => run.threads_per_process,
        :smts_per_thread => run.smts_per_thread,
    )
    for (key, val) in run.extra_run_params
        dict[Symbol(key)] = val
    end
    return dict
end
function Base.show(io::IO, run::RunParams)
    println(io, "RunParams:")
    println(io, "    machine: ", run.machine)
    println(io, "    queue: ", run.queue)
    println(io, "    walltime_seconds: ", run.walltime_seconds)
    println(io, "    nodes: ", run.nodes)
    println(io, "    cores_per_node: ", run.cores_per_node)
    println(io, "    pus_per_core: ", run.pus_per_core)
    println(io, "    processes: ", run.processes)
    println(io, "    threads_per_process: ", run.threads_per_process)
    println(io, "    smts_per_thread: ", run.smts_per_thread)
    println(io, "    extra_run_params:")
    for (key, val) in sort!(collect(run.extra_run_params))
        println(io, "        ", key, ": ", val)
    end
    return nothing
end
function Base.hash(run::RunParams, h::UInt)
    return hash(
        run.machine,
        hash(
            run.queue,
            hash(
                run.walltime_seconds,
                hash(
                    run.nodes,
                    hash(
                        run.cores_per_node,
                        hash(
                            run.pus_per_core,
                            hash(
                                run.processes,
                                hash(run.threads_per_process, hash(run.smts_per_thread, hash(run.extra_run_params, h))),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    )
end

export Benchmark
struct Benchmark
    physics::PhysicsParams
    run::RunParams
end
function Base.show(io::IO, benchmark::Benchmark)
    println(io, "Benchmark:")
    print(io, benchmark.physics)
    print(io, benchmark.run)
    return nothing
end
Base.hash(benchmark::Benchmark, h::UInt) = hash(benchmark.physics, hash(benchmark.run, h))

export Timing
struct Timing
    submitted::Bool
    success::Bool
    evolution_seconds::Float64
    evolution_compute_seconds::Float64
    evolution_output_seconds::Float64
    evolution_cell_updates::Float64
    evolution_iterations::Int
    cells::Float64
    cell_updates_per_second::Float64
    function Timing(;
        submitted::Bool,
        success::Bool,
        evolution_seconds::Float64=NaN,
        evolution_compute_seconds::Float64=NaN,
        evolution_output_seconds::Float64=NaN,
        evolution_cell_updates::Float64=NaN,
        evolution_iterations::Int=0,
    )
        if isnan(evolution_seconds) ||
            isnan(evolution_compute_seconds) ||
            isnan(evolution_output_seconds) ||
            isnan(evolution_cell_updates) ||
            evolution_iterations == 0
            @assert isnan(evolution_seconds) &&
                isnan(evolution_compute_seconds) &&
                isnan(evolution_output_seconds) &&
                isnan(evolution_cell_updates) &&
                evolution_iterations == 0
        end
        cells = evolution_cell_updates / evolution_iterations
        cell_updates_per_second = evolution_cell_updates / evolution_compute_seconds
        return new(
            submitted,
            success,
            evolution_seconds,
            evolution_compute_seconds,
            evolution_output_seconds,
            evolution_cell_updates,
            evolution_iterations,
            cells,
            cell_updates_per_second,
        )
    end
end
function Base.show(io::IO, timing::Timing)
    println(io, "Timing:")
    println(io, "    submitted: ", timing.submitted)
    println(io, "    success: ", timing.success)
    println(io, "    evolution_seconds: ", timing.evolution_seconds)
    println(io, "    evolution_compute_seconds: ", timing.evolution_compute_seconds)
    println(io, "    evolution_output_seconds: ", timing.evolution_output_seconds)
    println(io, "    evolution_cell_updates: ", timing.evolution_cell_updates)
    println(io, "    evolution_iterations: ", timing.evolution_iterations)
    println(io, "    number of cells (problem size): ", timing.cells)
    println(io, "    cell updates per second (strong performance): ", timing.cell_updates_per_second)
    return nothing
end
function Base.convert(::Type{DICT}, timing::Timing) where {DICT<:AbstractDict}
    return DICT(
        :submitted => timing.submitted,
        :success => timing.success,
        :evolution_seconds => timing.evolution_seconds,
        :evolution_compute_seconds => timing.evolution_compute_seconds,
        :evolution_output_seconds => timing.evolution_output_seconds,
        :evolution_cell_updates => timing.evolution_cell_updates,
        :evolution_iterations => timing.evolution_iterations,
        :cells => timing.cells,
        :cell_updates_per_second => timing.cell_updates_per_second,
    )
end

export BenchmarkResult
struct BenchmarkResult
    benchmark::Benchmark
    timing::Timing
end
function Base.show(io::IO, result::BenchmarkResult)
    println(io, "BenchmarkResult:")
    print(io, result.benchmark)
    print(io, result.timing)
    return nothing
end

################################################################################

export find_cactusdir
function find_cactusdir()
    dir = pwd()
    while true
        isdir(joinpath(dir, "simfactory")) && return dir
        @assert dir ≠ "/"
        dir = splitdir(dir)[1]
    end
    @assert false
end

function make_simulation_name(benchmark::Benchmark)
    parfile = benchmark.physics.parfile
    parfile = replace(parfile, r"^.*/" => s"")
    parfile = replace(parfile, r"(\.par)?$" => s"")

    settings = [
        [replace("$key$val", r"\$" => s"") for (key, val) in sort!(collect(benchmark.physics.extra_physics_params); by=first)]
        [replace("$key$val", r"\$" => s"") for (key, val) in sort!(collect(benchmark.run.extra_run_params); by=first)]
    ]
    settings::AbstractVector

    tag = bytes2hex(sha256(string(benchmark)))

    name1 = join(["benchmark", benchmark.physics.name, benchmark.physics.configuration, parfile, settings...], "-")
    name2 = join(
        [
            "",                 # empty string as placeholder for `name1`
            benchmark.run.machine,
            benchmark.run.queue,
            "n$(benchmark.run.nodes)",
            "c$(benchmark.run.cores_per_node)",
            "v$(benchmark.run.pus_per_core)",
            "p$(benchmark.run.processes)",
            "t$(benchmark.run.threads_per_process)",
            "s$(benchmark.run.smts_per_thread)",
            tag,
        ],
        "-",
    )

    # Limit length of simulation name
    maxlen = 200
    name1len = min(length(name1), maxlen - length(name2))
    simulation_name = name1[1:name1len] * name2
    @assert length(simulation_name) ≤ maxlen

    return simulation_name
end

export RunStatus, rs_unknown, rs_queued, rs_running, rs_finished
@enum RunStatus rs_unknown rs_queued rs_running rs_finished
function Base.show(io::IO, status::RunStatus)
    status ≡ rs_unknown && return print(io, "unknown")
    status ≡ rs_queued && return print(io, "queued")
    status ≡ rs_running && return print(io, "running")
    status ≡ rs_finished && return print(io, "finished")
    @assert false
    return nothing
end
Base.show(io::IO, ::MIME"text/plain", status::RunStatus) = show(io, status)

export get_run_status
function get_run_status(benchmark::Benchmark)
    cactusdir = find_cactusdir()
    name = make_simulation_name(benchmark)

    @info "Querying simulation $name..."
    output = read(
        Cmd(
            Cmd([joinpath("simfactory", "bin", "sim"), "--machine=$(benchmark.run.machine)", "list-simulations", name]);
            dir=cactusdir,
        ),
        String,
    )
    status = rs_unknown
    for line in split(output, "\n")
        if match(r"PRESUBMITTED|QUEUED", line) ≢ nothing
            status = rs_queued
            break
        elseif match(r"RUNNING", line) ≢ nothing
            status = rs_running
            break
        elseif match(r"FINISHED|INACTIVE", line) ≢ nothing
            status = rs_finished
            break
        end
    end
    @info "    $(string(status))"

    return status
end

struct NoCmd <: Base.AbstractCmd end
Base.show(io::IO, ::NoCmd) = print(io, "NoCmd()")
Base.:(==)(::NoCmd, ::NoCmd) = true
Base.hash(::NoCmd, h::UInt) = hash(0xdf3a71e1, h)
Base.ignorestatus(::NoCmd) = NoCmd()
Base.wait(::NoCmd) = nothing

function submit_run_nowait(benchmark::Benchmark)
    # Don't double-submit
    status = get_run_status(benchmark)
    status ≢ rs_unknown && return NoCmd()

    cactusdir = find_cactusdir()
    name = make_simulation_name(benchmark)

    replacements = [
        ["--replace=$key=$val" for (key, val) in sort!(collect(benchmark.physics.extra_physics_params); by=first)]
        ["--replace=$key=$val" for (key, val) in sort!(collect(benchmark.run.extra_run_params); by=first)]
    ]

    walltime_seconds = round(Int, benchmark.run.walltime_seconds)
    hours = walltime_seconds ÷ 3600
    minutes = walltime_seconds % 3600 ÷ 60
    seconds = walltime_seconds % 60
    walltime = "$hours:$minutes:$seconds"

    nodes = benchmark.run.nodes
    cores_per_node = benchmark.run.cores_per_node
    pus_per_core = benchmark.run.pus_per_core
    processes = benchmark.run.processes
    threads_per_process = benchmark.run.threads_per_process
    smts_per_thread = benchmark.run.smts_per_thread

    threads = processes * threads_per_process * smts_per_thread
    threads_per_node = threads ÷ nodes
    cores = nodes * cores_per_node
    threads_per_core = threads ÷ cores

    @info "Submitting $name..."
    cmd = run(
        Cmd(
            Cmd([
                joinpath("simfactory", "bin", "sim"),
                "--machine=$(benchmark.run.machine)",
                "submit",
                name,
                "--configuration=$(benchmark.physics.configuration)",
                "--parfile=$(benchmark.physics.parfile)",
                replacements...,
                "--queue=$(benchmark.run.queue)",
                "--walltime=$walltime",
                "--ppn=$cores_per_node",
                "--procs=$threads",
                "--ppn-used=$threads_per_node",
                "--num-threads=$threads_per_process",
                "--num-smt=$threads_per_core",
            ]);
            dir=cactusdir,
        );
        wait=false,
    )

    return cmd
end

export submit_run
function submit_run(benchmark::Benchmark)
    cmd = submit_run_nowait(benchmark::Benchmark)
    @info "Waiting for all submission to finish..."
    wait(cmd)
    @info "Done."
    return nothing
end

export submit_runs
function submit_runs(benchmarks::AbstractVector{Benchmark})
    cmds = submit_run_nowait.(benchmarks)
    @info "Waiting for all submissions to finish..."
    wait.(cmds)
    @info "Done."
    return nothing
end

export wait_for_run
function wait_for_run(benchmark::Benchmark)
    name = make_simulation_name(benchmark)
    @info "Waiting for $name..."
    delay_seconds = 1
    max_delay_seconds = 60
    while true
        status = get_run_status(benchmark)
        status ≡ rs_unknown && return status
        status ≡ rs_finished && return status
        @info "    Waiting $delay_seconds seconds..."
        sleep(delay_seconds)
        # Back off exponentially
        delay_seconds = min(max_delay_seconds, 2 * delay_seconds)
    end
end

export wait_for_runs
function wait_for_runs(benchmarks::AbstractVector{Benchmark})
    for benchmark in benchmarks
        status = wait_for_run(benchmark)
        @assert status ≡ rs_finished
    end
    return nothing
end

export read_run_timing
function read_run_timing(benchmark::Benchmark)
    status = wait_for_run(benchmark)
    if status ≡ rs_unknown
        # Run does not exist
        return Timing(; submitted=false, success=false)
    end

    cactusdir = find_cactusdir()
    name = make_simulation_name(benchmark)

    cmd = Cmd(
        Cmd([joinpath("simfactory", "bin", "sim"), "--machine=$(benchmark.run.machine)", "get-output-dir", name]); dir=cactusdir
    )
    output_dir = chomp(read(cmd, String))

    parfile = benchmark.physics.parfile
    parfile = replace(parfile, r"^.*/" => s"")
    parfile = replace(parfile, r"(\.par)?$" => s"")

    timing = try
        filename = joinpath(output_dir, parfile, "performance.yaml")
        yaml = YAML.load_file(filename)
        performance = yaml["performance"]::Dict
        Timing(;
            submitted=true,
            success=true,
            evolution_seconds=performance["evolution-seconds"],
            evolution_compute_seconds=performance["evolution-compute-seconds"],
            evolution_output_seconds=performance["evolution-output-seconds"],
            evolution_cell_updates=performance["evolution-cell-updates"],
            evolution_iterations=performance["evolution-iterations"],
        )
    catch
        # Could not read output; run probably failed
        return Timing(; submitted=true, success=false)
    end

    return timing
end

end
