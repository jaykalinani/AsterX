FROM einsteintoolkit/cactus-amrex
WORKDIR /

RUN useradd -m jovyan
USER jovyan
WORKDIR /home/jovyan
ENV HOME /home/jovyan
ENV USER jovyan
RUN mkdir -p /home/jovyan/bin
USER root
COPY build.sh /usr/local/bin/
RUN chmod 755 /usr/local/bin/build.sh
#ENV container docker
COPY Dockerfile /Dockerfile
USER jovyan
CMD ["sleep","infinity"]
