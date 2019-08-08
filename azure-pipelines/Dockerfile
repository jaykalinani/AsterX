# Use Ubuntu
FROM ubuntu:19.04

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
RUN apt-get update && \
    apt-get -y --no-install-recommends install apt-utils build-essential ca-certificates g++ git gfortran libopenmpi-dev perl pkg-config python rsync subversion wget && \
    rm -rf /var/lib/apt/lists/*

# Install cmake
# (AMReX 19.07 requires at least cmake 3.14)
RUN wget https://github.com/Kitware/CMake/releases/download/v3.14.5/cmake-3.14.5-Linux-x86_64.tar.gz && \
    tar xzf cmake-3.14.5-Linux-x86_64.tar.gz -C /usr/local --strip-components=1 && \
    rm cmake-3.14.5-Linux-x86_64.tar.gz

# Build AMReX
RUN wget https://github.com/AMReX-Codes/amrex/archive/19.08.tar.gz && \
    tar xzf 19.08.tar.gz && \
    cd amrex-19.08 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_ASSERTIONS=ON -DENABLE_BACKTRACE=ON -DENABLE_OMP=ON -DCMAKE_INSTALL_PREFIX=/cactus/amrex .. && \
    make -j2 && \
    make -j2 install && \
    cd ../.. && \
    rm -rf 19.08.tar.gz amrex-19.08

# As documentation
COPY Dockerfile /Dockerfile
