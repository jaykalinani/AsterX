# Use Ubuntu
FROM ubuntu:19.04

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
RUN apt-get update && apt-get -y install build-essential g++ git gfortran libopenmpi-dev perl pkg-config python subversion wget && rm -rf /var/lib/apt/lists/*

# Install cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.14.5/cmake-3.14.5-Linux-x86_64.tar.gz && tar xzf cmake-3.14.5-Linux-x86_64.tar.gz -C /usr/local --strip-components=1 && rm cmake-3.14.5-Linux-x86_64.tar.gz

# Build AMReX
RUN wget https://github.com/AMReX-Codes/amrex/archive/19.07.tar.gz && tar xzf 19.07.tar.gz && cd amrex-19.07 && mkdir build && cd build && cmake -DENABLE_OMP=ON -DCMAKE_INSTALL_PREFIX=/cactus/amrex .. && make -j2 && make -j2 install && cd ../.. && rm -rf amrex-19.07

# As documentation
COPY Dockerfile /Dockerfile
