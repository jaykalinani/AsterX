#!/bin/bash

################################################################################

{
    iter=64
    # 256 320 400 512
    for size in 256; do
        for nodes in 1; do          # total
            cores=40                # per node
            procs=$((nodes*cores))  # total
            # 1 2 4 5 10 20 40
            # for threads in 5 10 20; do # per process
            for threads in 5; do # per process
    
                for blocking_factor in 8; do
                    # 32 64
                    for max_grid_size in 64; do
                        # 4 8 16 32 64 1024000
                        # for max_tile_size_x in 16 1024000; do
                        for max_tile_size_x in 1024000; do
                            # 4 8 16 32 64
                            # for max_tile_size_y in 2 4 8; do
                            for max_tile_size_y in 4; do
                                # for max_tile_size_z in 2 4 8; do
                                for max_tile_size_z in 2; do
    
                                    name="benchmark-minkowski-unigrid-i$iter-s$size-n$nodes-c$cores-p$procs-t$threads-bf$blocking_factor-gs$max_grid_size-tx$max_tile_size_x-ty$max_tile_size_y-tz$max_tile_size_z"
                                    ./simfactory/bin/sim purge "$name"
                                    ./simfactory/bin/sim \
                                        submit "$name" \
                                        --parfile arrangements/CarpetX/Z4c/par/benchmark-minkowski-unigrid.par \
                                        --replace='$ncells='$size \
                                        --replace='$blocking_factor='$blocking_factor \
                                        --replace='$max_grid_size='$max_grid_size \
                                        --replace='$max_tile_size_x='$max_tile_size_x \
                                        --replace='$max_tile_size_y='$max_tile_size_y \
                                        --replace='$max_tile_size_z='$max_tile_size_z \
                                        --procs=$procs \
                                        --num-threads=$threads \
                                        --queue=debugq \
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
    echo 'system,iter,time,size,nodes,cores_node,threads,threads_process,blocking_factor,max_grid_size,max_tile_size_x,max_tile_size_y,max_tile_size_z'
    system=$(./simfactory/bin/sim whoami | awk '{ print $3; }')
    (
        cd /home/eschnetter/simulations;
        for file in benchmark-minkowski-unigrid-i*/output-0000/benchmark-minkowski-unigrid-*.out; do
            echo \
                "$system" \
                $(grep Evolve "$file" |
                      tail -n 1 |
                      awk '{ print $2; }') \
                          $(basename "$file" |
                                sed -e 's/benchmark-minkowski-unigrid//' |
                                sed -e 's/\.out//')
        done
    ) |
        sort -k 2 -n |
        sed -e 's/ -i/,/;s/-s/,/;s/-n/,/;s/-c/,/;s/-p/,/;s/-t/,/;s/-bf/,/;s/-gs/,/;s/-tx/,/;s/-ty/,/;s/-tz/,/'
} |
    tee benchmark-minkowski-unigrid.csv

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
