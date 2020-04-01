# How to build this Docker image:
#     docker build . --tag eschnett/cactus-amrex
#     docker push eschnett/cactus-amrex

# Use Ubuntu
FROM ubuntu:19.10

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
RUN apt-get update && \
    apt-get --yes --no-install-recommends install \
        build-essential \
        ca-certificates \
        g++ \
        gfortran \
        git \
        hdf5-tools \
        libfftw3-dev \
        libgsl-dev \
        libhdf5-dev \
        libhwloc-dev \
        libmpich-dev \
        libopenblas-dev \
        perl \
        pkg-config \
        python \
        python3 \
        rsync \
        subversion \
        wget \
        zlib1g-dev \
        && \
    rm -rf /var/lib/apt/lists/*

# Install cmake
# (AMReX 20.04 requires at least cmake 3.14)
RUN mkdir dist && \
    (cd dist && \
    wget https://github.com/Kitware/CMake/releases/download/v3.16.2/cmake-3.16.2-Linux-x86_64.tar.gz && \
    tar xzf cmake-3.16.2-Linux-x86_64.tar.gz -C /usr/local --strip-components=1 && \
    true) && \
    rm -rf dist

# Install NSIMD
# Note: This assumes that the system has x86_64 CPUs with AVX2 and FMA
RUN mkdir src && \
    (cd src && \
    git clone -b eschnett/storeu_masked https://github.com/eschnett/nsimd.git && \
    cd nsimd && \
    python3 egg/hatch.py --all --force --disable-clang-format && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DSIMD=AVX2 -DSIMD_OPTIONALS=FMA -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install ssht
RUN mkdir src && \
    (cd src && \
    git clone -n https://github.com/astro-informatics/ssht.git && \
    cd ssht && \
    git checkout 7378ce8853897cbd1b08adebf7ec088c1e40f860 && \
    make -j$(nproc) default && \
    cp src/c/*.h /usr/local/include && \
    cp lib/c/* /usr/local/lib && \
    true) && \
    rm -rf src

# Install Silo
RUN mkdir src && \
    (cd src && \
    wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz && \
    tar xzf silo-4.10.2-bsd.tar.gz && \
    cd silo-4.10.2-bsd && \
    mkdir build && \
    cd build && \
    ../configure --disable-fortran --enable-optimization --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial/include,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib --prefix=/usr/local && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install AMReX
RUN mkdir src && \
    (cd src && \
    wget https://github.com/AMReX-Codes/amrex/archive/20.04.tar.gz && \
    tar xzf 20.04.tar.gz && \
    cd amrex-20.04 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_PARTICLES=ON -DENABLE_ASSERTIONS=ON -DENABLE_BACKTRACE=ON -DENABLE_OMP=ON -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# As documentation
COPY Dockerfile /Dockerfile
