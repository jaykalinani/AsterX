# Creating Docker images for CarpetX (with CPUs)

## Create a Docker image containing Spack
You don't need to do this; these commands are for maintainers.
```sh
cd $HOME/src/spack
docker build --file share/spack/docker/ubuntu-2004.dockerfile --build-arg baseimage=ubuntu:20.04 --tag eschnett/spack-ubuntu-2004:cpu .
# docker save -o spack-ubuntu-2004-cpu.tar eschnett/spack-ubuntu-2004:cpu
docker push eschnett/spack-ubuntu-2004:cpu
```

## Create a Docker image containing AMReX and other Cactus dependencies
You don't need to do this; these commands are for maintainers.
```sh
cd utils/Docker/make-image-cpu
docker build --tag eschnett/simfactory4-cactus:cpu .
docker push eschnett/simfactory4-cactus:cpu
```

## Use the Docker image to build and run Cactus
This is what you need to do.
```sh
./utils/Docker/make/image-cpu/run-container-cpu
```
This opens a shell running in a container using the image above.
