# Build stage with Spack pre-installed and ready to be used

FROM eschnett/spack-ubuntu-2004:cpu as builder
MAINTAINER Erik Schnetter <schnetter@gmail.com>

# Install system packages
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get --yes --no-install-recommends install \
        bzip2 \
        cmake \
        libboost-all-dev \
        libbz2-dev \
        libfftw3-dev \
        libgsl-dev \
        libhwloc-dev \
        libjemalloc-dev \
        libjpeg-dev \
        libncurses-dev \        
        libopenblas-dev \
        libopenmpi-dev \
        libpciaccess-dev \
        libsilo-dev \
        libssl-dev \
        libxml2-dev \
        libyaml-cpp-dev \
        numactl \
        perl \
        pkg-config \
        python3 \
        swig \
        zlib1g-dev \
        && \
    rm -rf /var/lib/apt/lists/*

COPY packages.yaml /opt/spack/etc/spack/packages.yaml

# What we want to install and how we want to install it is specified
# in a manifest file (spack.yaml)
COPY spack.yaml /opt/spack-environment/spack.yaml

# Install the software, remove unecessary deps
RUN cd /opt/spack-environment && \
    spack env activate . && \
    spack install --fail-fast && \
    spack gc -y

# Strip all the binaries (this saves a lot of space)
RUN find -L /opt/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    tee >(grep 'x-executable\|x-sharedlib' | \
          awk -F: '{print $1}' | xargs strip --strip-all) \
        >(grep 'x-archive' | \
          awk -F: '{print $1}' | xargs strip) \
        > /dev/null

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh



# Build on the standard Ubuntu image
FROM ubuntu:20.04

# Install system packages
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get --yes --no-install-recommends install \
        build-essential \
        bzip2 \
        ca-certificates \
        file \
        g++ \
        gcc \
        gfortran \
        git \
        hdf5-tools \
        libboost-all-dev \
        libbz2-dev \
        libfftw3-dev \
        libgsl-dev \
        libhwloc-dev \
        libjemalloc-dev \
        libjpeg-dev \
        libopenblas-dev \
        libopenmpi-dev \
        libpciaccess-dev \
        libsilo-dev \
        libssl-dev \
        libxml2-dev \
        libyaml-cpp-dev \
        numactl \
        perl \
        pkg-config \
        python \
        python3 \
        rsync \
        slurm-client \
        ssh \
        subversion \
        wget \
        zlib1g-dev \
        && \
    rm -rf /var/lib/apt/lists/*

# Install Spack packages
COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/view /opt/view
COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh
