# Use Ubuntu
FROM fedora

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
RUN dnf install -y curl wget gcc-g++ gfortran make git vim
RUN dnf install -y openmpi openmpi-devel

ENV PATH /usr/lib64/openmpi/bin:/usr/local/bin:/usr/bin:/usr/sbin:/bin:/sbin:.

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
#COPY Dockerfile /Dockerfile

RUN dnf install -y svn perl-Thread-Queue
RUN dnf install -y python2
RUN curl -kLO https://raw.githubusercontent.com/nds-org/jupyter-et/master/et-pkg-installer.py
RUN python2 et-pkg-installer.py
RUN sh install-for-cactus.sh
RUN dnf install -y file
RUN useradd -m jovyan
USER jovyan
WORKDIR /home/jovyan
ENV HOME /home/jovyan
ENV USER jovyan
RUN mkdir -p /home/jovyan/bin
USER root
COPY amrex.cfg /amrex.cfg
COPY setmeup.sh /usr/local/bin/
RUN chmod 755 /usr/local/bin/setmeup.sh
RUN dnf install -y gdb
#ENV container docker
COPY Dockerfile /Dockerfile
RUN dnf install -y diffutils
USER jovyan
#CMD ["/usr/sbin/init"]
CMD ["sleep","infinity"]
