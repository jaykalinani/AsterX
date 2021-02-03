#!/bin/bash

################################################################################

# 256 320 400 512
for size in 256; do
    for nodes in 4; do                  # total
        cores=40                              # per node
        procs=$((nodes*cores))                # total
        # for threads in 1 2 4 5 8 10 20 40; do # per process
        for threads in 40; do # per process
            name="benchmark-minkowski-unigrid-s$size-n$nodes-c$cores-p$procs-t$threads"
            ./simfactory/bin/sim purge "$name"
            ./simfactory/bin/sim submit "$name" --parfile arrangements/CarpetX/Z4c/par/benchmark-minkowski-unigrid.par --replace='$ncells='"$size" --procs=$procs --num-threads=$threads --walltime=1:0:0 &
        done
    done
done
wait

################################################################################

# Analyse output
{
    echo 'timesize,nodes,cores_node,threads,threads_process'
    (
        cd /home/eschnetter/simulations;
        for file in benchmark-minkowski-unigrid-*/output-0000/benchmark-minkowski-unigrid-*.out; do
            echo \
                $(grep Evolve "$file" |
                      tail -n +3 |
                      awk '{ print $2; }') \
                          $(basename "$file" |
                                sed -e 's/benchmark-minkowski-unigrid-/s256-/' |
                                sed -e 's/s256-s/s/' |
                                sed -e 's/\.out//')
        done
    ) |
        sort -n |
        sed -e 's/ s/,/;s/-n/,/;s/-c/,/;s/-p/,/;s/-t/,/'
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
