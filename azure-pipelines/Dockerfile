# How to build this Docker image:
#     docker build . --tag einsteintoolkit/cactus-amrex
#     docker push einsteintoolkit/cactus-amrex

# Use Ubuntu
# FROM ubuntu:20.04
FROM debian:testing-20200607-slim

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
# Note: Boost on Ubuntu requires OpenMPI
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
        libudev-dev \
        perl \
        pkg-config \
        python \
        python3 \
        rsync \
        subversion \
        vim \
        wget \
        zlib1g-dev \
        && \
    rm -rf /var/lib/apt/lists/*

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

# Install ssht
# TODO: Use "cmake install" instead of "cp"
RUN mkdir src && \
    (git clone https://github.com/astro-informatics/ssht.git && \
    mkdir ssht/build && \
    cd ssht/build && \
    cmake .. && \
    make install) && \
    rm -fr src

# Install yaml-cpp
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
RUN mkdir src && \
    (cd src && \
    wget https://github.com/AMReX-Codes/amrex/archive/20.06.tar.gz && \
    tar xzf 20.06.tar.gz && \
    cd amrex-20.06 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_PARTICLES=ON -DENABLE_ASSERTIONS=ON -DENABLE_BACKTRACE=ON -DENABLE_OMP=ON -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# # Install AMReX
# RUN mkdir src && \
#     (cd src && \
#     git clone -n https://github.com/mic84/amrex && \
#     cd amrex && \
#     git checkout c01d9965064b20972e2109e85f842ecb1d712c41 && \
#     mkdir build && \
#     cd build && \
#     cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_PARTICLES=ON -DENABLE_ASSERTIONS=ON -DENABLE_BACKTRACE=ON -DENABLE_HDF5=ON -DENABLE_OMP=ON -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
#     make -j$(nproc) && \
#     make -j$(nproc) install && \
#     true) && \
#     rm -rf src
