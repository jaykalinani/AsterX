# How to build this Docker image:
#     docker build . --tag eschnett/cactus-amrex
#     docker push eschnett/cactus-amrex

# Use Ubuntu
FROM ubuntu:19.04

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
        libopenmpi-dev \
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
# (AMReX 19.12 requires at least cmake 3.14)
RUN mkdir dist && \
    ( cd dist && \
    wget https://github.com/Kitware/CMake/releases/download/v3.14.5/cmake-3.14.5-Linux-x86_64.tar.gz && \
    tar xzf cmake-3.14.5-Linux-x86_64.tar.gz -C /usr/local --strip-components=1 && \
    true ) && \
    rm -rf dist

# Install NSIMD
# Note: This assumes that the system has x86_64 CPUs with AVX2 and FMA
RUN mkdir src && \
    ( cd src && \
    git clone -b eschnett/storeu_masked https://github.com/eschnett/nsimd.git && \
    cd nsimd && \
    python3 egg/hatch.py --all --force --disable-clang-format && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DSIMD=AVX2 -DSIMD_OPTIONALS=FMA -DCMAKE_INSTALL_PREFIX=/cactus/nsimd .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true ) && \
    rm -rf src

# Install AMReX
RUN mkdir src && \
    ( cd src && \
    wget https://github.com/AMReX-Codes/amrex/archive/19.12.tar.gz && \
    tar xzf 19.12.tar.gz && \
    cd amrex-19.12 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_PARTICLES=ON -DENABLE_ASSERTIONS=ON -DENABLE_BACKTRACE=ON -DENABLE_OMP=ON -DCMAKE_INSTALL_PREFIX=/cactus/amrex .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true ) && \
    rm -rf src

# As documentation
COPY Dockerfile /Dockerfile
