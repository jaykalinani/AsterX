# How to build this Docker image:
#     docker build . --tag einsteintoolkit/cactus-amrex
#     docker push einsteintoolkit/cactus-amrex

# Use Debian
# - We need a rather recent Debian to get a recent cmake for AMReX
FROM debian:testing-20210621-slim

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
# - Boost on Ubuntu requires OpenMPI
# TODO: Add CVS to check out Lorene
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get --yes --no-install-recommends install \
        build-essential \
        ca-certificates \
        cmake \
        g++ \
        gfortran \
        git \
        hdf5-tools \
        libboost-all-dev \
        libfftw3-dev \
        libgsl-dev \
        libhdf5-dev \
        libhwloc-dev \
        libopenblas-dev \
        libopenmpi-dev \
        libpetsc-real-dev \
        libudev-dev \
        perl \
        pkg-config \
        python2 \
        python3 \
        rsync \
        subversion \
        vim \
        wget \
        zlib1g-dev \
        && \
    rm -rf /var/lib/apt/lists/*

# Install ADIOS2
# ADIOS2 is a parallel I/O library, comparable to HDF5
RUN mkdir src && \
    (cd src && \
    wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.7.1.tar.gz && \
    tar xzf v2.7.1.tar.gz && \
    cd ADIOS2-2.7.1 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install NSIMD
# NSIMD allows writing explicitly SIMD-vectorized code
# Note: This assumes that the system has x86_64 CPUs with AVX2
RUN mkdir src && \
    (cd src && \
    wget https://github.com/agenium-scale/nsimd/archive/refs/tags/v2.2.tar.gz && \
    tar xzf v2.2.tar.gz && \
    cd nsimd-2.2 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -Dsimd=AVX2 -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install openPMD-api
# openPMD-api defines a standard for laying out AMR data in a file
# - Depends on ADIOS2
RUN mkdir src && \
    (cd src && \
    wget https://github.com/openPMD/openPMD-api/archive/refs/tags/0.14.0.tar.gz && \
    tar xzf 0.14.0.tar.gz && \
    cd openPMD-api-0.14.0 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install Silo
# openPMD-api defines a standard for laying out AMR data in a file
RUN mkdir src && \
    (cd src && \
    wget https://wci.llnl.gov/sites/wci/files/2021-01/silo-4.10.2-bsd.tgz && \
    tar xzf silo-4.10.2-bsd.tgz && \
    cd silo-4.10.2-bsd && \
    mkdir build && \
    cd build && \
    ../configure --disable-fortran --enable-optimization --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial/include,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib --prefix=/usr/local && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install ssht
# ssht provides spin-weighted spherical harmonicsa
# TODO: Use "cmake install" instead of "cp"
# RUN mkdir src && \
#     (git clone https://github.com/astro-informatics/ssht.git && \
#     mkdir ssht/build && \
#     cd ssht/build && \
#     cmake .. && \
#     make install) && \
#     rm -fr src
RUN mkdir src && \
    (cd src && \
    wget https://github.com/astro-informatics/ssht/archive/v1.5.0.tar.gz && \
    tar xzf v1.5.0.tar.gz && \
    cd ssht-1.5.0 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install yaml-cpp
# yaml-cpp allows reading and writing YAML files
RUN mkdir src && \
    (cd src && \
    wget https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.3.tar.gz && \
    tar xzf yaml-cpp-0.6.3.tar.gz && \
    cd yaml-cpp-yaml-cpp-0.6.3 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install AMReX
# AMReX provides adaptive mesh refinement
# - Enable Fortran for `docker/Dockerfile`
# - Keep AMReX source tree around for debugging
# - Install this last because it changes most often
RUN mkdir src && \
    (cd src && \
    wget https://github.com/AMReX-Codes/amrex/archive/21.08.tar.gz && \
    tar xzf 21.08.tar.gz && \
    cd amrex-21.08 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Debug -DAMReX_OMP=ON -DAMReX_PARTICLES=ON -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true)
