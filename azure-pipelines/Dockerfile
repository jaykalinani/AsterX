# Use Ubuntu
FROM ubuntu:19.04

# Add user "cactus" with sudo privileges
RUN apt-get update && apt-get -y install sudo && rm -rf /var/lib/apt/lists/*
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers
# Don't use /home since this can lead to ambiguities with Singularity
RUN useradd -G sudo -m -s /bin/bash -d /cactus cactus
USER cactus
WORKDIR /cactus

# Install system packages
RUN sudo apt-get update && sudo apt-get -y install build-essential g++ gfortran  libopenmpi-dev perl pkg-config python wget && sudo rm -rf /var/lib/apt/lists/*

# Install cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.14.5/cmake-3.14.5-Linux-x86_64.tar.gz && sudo tar xzf cmake-3.14.5-Linux-x86_64.tar.gz -C /usr/local --strip-components=1

# Build AMReX
RUN wget https://github.com/AMReX-Codes/amrex/archive/19.07.tar.gz && tar xzf 19.07.tar.gz && cd amrex-19.07 && mkdir build && cd build && cmake -DENABLE_OMP=ON -DCMAKE_INSTALL_PREFIX=/cactus/amrex .. && make -j2 && make -j2 install

# Set up Simfactory
RUN echo 'docker' >$HOME/.hostname && mkdir /cactus/simulations
COPY defs.local.ini /cactus/defs.local.ini
COPY docker.ini /cactus/docker.ini
COPY ubuntu.cfg /cactus/ubuntu.cfg
COPY docker.sub /cactus/docker.sub
COPY docker.run /cactus/docker.run
COPY amrex.th /cactus/amrex.th

# As documentation
COPY Dockerfile /Dockerfile

# Create user
RUN sudo useradd -u 1001 -G sudo -m -s /bin/bash -d /vsts vsts
