# How to run this Docker image:
#     docker build . --build-arg commit=<commit>
#     docker build . --build-arg commit=$(git log -n 1 | awk '/^commit / { print $2; }')

FROM eschnett/cactus-amrex

# Download repo
ARG commit
RUN git clone -n https://bitbucket.org/eschnett/cactusamrex.git && \
    cd cactusamrex && \
    git checkout "${commit}"
WORKDIR cactusamrex

# Download Cactus
RUN wget https://raw.githubusercontent.com/gridaphobe/CRL/ET_2019_10/GetComponents && \
    chmod a+x GetComponents && \
    ./GetComponents --parallel azure-pipelines/carpetx.th && \
    echo 'docker' >$HOME/.hostname && \
    cp azure-pipelines/defs.local.ini Cactus/simfactory/etc && \
    cp azure-pipelines/docker.ini Cactus/simfactory/mdb/machines && \
    cp azure-pipelines/ubuntu.cfg Cactus/simfactory/mdb/optionlists && \
    cp azure-pipelines/docker.sub Cactus/simfactory/mdb/submitscripts && \
    cp azure-pipelines/docker.run Cactus/simfactory/mdb/runscripts && \
    mkdir $HOME/simulations
WORKDIR Cactus

# Build Cactus
RUN ./simfactory/bin/sim build -j$(nproc) --thornlist=../azure-pipelines/carpetx.th

# Run example
RUN env OMP_NUM_THREADS=2 ./simfactory/bin/sim submit ./simfactory/etc/parfiles/empty.par && \
    ./simfactory/bin/sim show-output empty && \
    env OMP_NUM_THREADS=2 mpiexec -n 2 ./exe/cactus_sim ../azure-pipelines/carpetx.par

# Run tests
RUN ./simfactory/bin/sim submit sim-test-1 --test --procs=1 --num-threads=2 && \
    cat /home/vsts_azpcontainer/simulations/sim-test-1/output-0000/sim-test-1.log && \
    ./simfactory/bin/sim submit sim-test-2 --test --procs=2 --num-threads=1 && \
    cat /home/vsts_azpcontainer/simulations/sim-test-2/output-0000/sim-test-2.log
