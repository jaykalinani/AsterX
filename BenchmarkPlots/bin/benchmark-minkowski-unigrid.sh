#!/bin/bash

################################################################################

{
    machine=symmetry
    configuration=sim
    queue=defq
    iter=9999
    # machine=symmetry-llvm
    # configuration=sim-llvm
    # queue=defq
    # iter=0001
    # 256 320 400 512
    # machine=symmetry-gpu
    # configuration=sim-gpu
    # queue=gpudebugq
    # iter=0002
    # 256 320 400 512
    for size in 256; do
        # 1 2 4 8 16 32
        for nodes in 1 2 4 8; do  # total
            cores=40            # per node
            # cores=16                # per node
            procs=$((nodes*cores)) # total
            # 1 2 4 5 10 20 40
            for threads in 10; do # per process

                for blocking_factor in 8; do
                    # 32 64 128
                    for max_grid_size in 128; do
                        # 4 8 16 32 64 1024000
                        # for max_tile_size_x in 16 1024000; do
                        for max_tile_size_x in 1024000; do
                            # 4 8 16 32 64
                            # for max_tile_size_y in 1 2 4 8 16; do
                            for max_tile_size_y in 4; do
                                # for max_tile_size_z in 1 2 4 8 16 32; do
                                for max_tile_size_z in 2; do

                                    name="benchmark-minkowski-unigrid-i$iter-s$size-n$nodes-c$cores-p$procs-t$threads-bf$blocking_factor-gs$max_grid_size-tx$max_tile_size_x-ty$max_tile_size_y-tz$max_tile_size_z"
                                    ./simfactory/bin/sim purge "$name"
                                    ./simfactory/bin/sim \
                                        --machine=$machine \
                                        submit "$name" \
                                        --configuration=$configuration \
                                        --parfile arrangements/CarpetX/Z4c/par/BenchmarkPlots/par/benchmark-minkowski-unigrid.par \
                                        --replace='$ncells='$size \
                                        --replace='$blocking_factor='$blocking_factor \
                                        --replace='$max_grid_size='$max_grid_size \
                                        --replace='$max_tile_size_x='$max_tile_size_x \
                                        --replace='$max_tile_size_y='$max_tile_size_y \
                                        --replace='$max_tile_size_z='$max_tile_size_z \
                                        --procs=$procs \
                                        --num-threads=$threads \
                                        --queue=$queue \
                                        --walltime=1:0:0 &

                                done
                            done
                        done
                    done
                done

            done
        done
    done
    wait
}

{
    machine=db
    configuration=sim
    iter=0000
    # 256 320 400 512
    for size in 256; do
        # 1 2 4 8 16 32
        for nodes in 1 2 4; do     # total
            cores=48               # per node
            procs=$((nodes*cores)) # total
            for threads in 12; do # per process

                for blocking_factor in 8; do
                    # 32 64 128
                    for max_grid_size in 128; do
                        # 4 8 16 32 64 1024000
                        # for max_tile_size_x in 16 1024000; do
                        for max_tile_size_x in 1024000; do
                            # 4 8 16 32 64
                            # for max_tile_size_y in 1 2 4 8 16; do
                            for max_tile_size_y in 16; do
                                # for max_tile_size_z in 1 2 4 8 16 32; do
                                for max_tile_size_z in 32; do

                                    name="benchmark-minkowski-unigrid-i$iter-s$size-n$nodes-c$cores-p$procs-t$threads-bf$blocking_factor-gs$max_grid_size-tx$max_tile_size_x-ty$max_tile_size_y-tz$max_tile_size_z"
                                    ./simfactory/bin/sim purge "$name"
                                    ./simfactory/bin/sim \
                                        --machine=$machine \
                                        submit "$name" \
                                        --configuration=$configuration \
                                        --parfile arrangements/CarpetX/Z4c/par/benchmark-minkowski-unigrid.par \
                                        --replace='$ncells='$size \
                                        --replace='$blocking_factor='$blocking_factor \
                                        --replace='$max_grid_size='$max_grid_size \
                                        --replace='$max_tile_size_x='$max_tile_size_x \
                                        --replace='$max_tile_size_y='$max_tile_size_y \
                                        --replace='$max_tile_size_z='$max_tile_size_z \
                                        --procs=$procs \
                                        --num-threads=$threads \
                                        --walltime=1:0:0 &

                                done
                            done
                        done
                    done
                done

            done
        done
    done
    wait
}

{
    machine=summit-gpu
    configuration=sim-gpu
    iter=0002
    # 256 320 400 512
    for size in 256 512; do
        # 1 2 4 8 16 32
        for nodes in 8 16 32; do         # total
            cores=36                   # per node
            smt=4                      # per core
            procs=$((nodes*cores*smt)) # total
            for threads in 24; do      # per process
    
                for blocking_factor in 8; do
                    # 32 64 128
                    for max_grid_size in 128; do
                        # 4 8 16 32 64 1024000
                        # for max_tile_size_x in 16 1024000; do
                        for max_tile_size_x in 1024000; do
                            # 4 8 16 32 64
                            # for max_tile_size_y in 1 2 4 8 16; do
                            for max_tile_size_y in 16; do
                                # for max_tile_size_z in 1 2 4 8 16 32; do
                                for max_tile_size_z in 32; do

                                    name="benchmark-minkowski-unigrid-i$iter-s$size-n$nodes-c$cores-p$procs-t$threads-bf$blocking_factor-gs$max_grid_size-tx$max_tile_size_x-ty$max_tile_size_y-tz$max_tile_size_z"
                                    ./simfactory/bin/sim purge "$name"
                                    ./simfactory/bin/sim \
                                        --machine=$machine \
                                        submit "$name" \
                                        --configuration=$configuration \
                                        --parfile arrangements/CarpetX/Z4c/par/benchmark-minkowski-unigrid.par \
                                        --replace='$ncells='$size \
                                        --replace='$blocking_factor='$blocking_factor \
                                        --replace='$max_grid_size='$max_grid_size \
                                        --replace='$max_tile_size_x='$max_tile_size_x \
                                        --replace='$max_tile_size_y='$max_tile_size_y \
                                        --replace='$max_tile_size_z='$max_tile_size_z \
                                        --ppn-used=$((cores*smt)) \
                                        --procs=$procs \
                                        --num-threads=$threads \
                                        --num-smt=$smt \
                                        --walltime=1:0:0 &

                                done
                            done
                        done
                    done
                done
    
            done
        done
    done
    wait
}

################################################################################

# Analyse output
{
    echo 'system,time,iter,size,nodes,cores_node,threads,threads_process,blocking_factor,max_grid_size,max_tile_size_x,max_tile_size_y,max_tile_size_z'
    system=$(./simfactory/bin/sim whoami | awk '{ print $3; }')
    (
        cd $HOME/simulations;
        for dir in benchmark-minkowski-unigrid-i*; do
            time=$(grep Evolve "$dir/output-0000/stdout.txt" |
                      tail -n 1 |
                      awk '{ print $2; }')
            settings=$(echo "$dir" |
                           sed -e 's/benchmark-minkowski-unigrid//')
            echo "$system,$time$settings"
        done
    ) |
        sort -k2 -n -t, |
        sed -e 's/-i/,/;s/-s/,/;s/-n/,/;s/-c/,/;s/-p/,/;s/-t/,/;s/-bf/,/;s/-gs/,/;s/-tx/,/;s/-ty/,/;s/-tz/,/'
} |
    tee benchmark-minkowski-unigrid.csv

################################################################################

# Look for errors
(
    cd /home/eschnetter/simulations;
    for file in benchmark-minkowski-unigrid-*/output-0000/benchmark-minkowski-unigrid-*.out; do
        echo \
            $(grep Evolve "$file" |
                  tail -n +3 |
                  wc -l) \
                      $(basename "$file" |
                            sed -e 's/benchmark-minkowski-unigrid-//' |
                            sed -e 's/\.out//')
    done
) |
    grep -v '^1'
