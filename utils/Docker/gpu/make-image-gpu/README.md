# Creating and using Docker images for CarpetX (with GPUs)

## Create a Docker image containing Spack
You don't need to do this; these commands are for maintainers.
```sh
cd $HOME/src/spack
docker build --file share/spack/docker/ubuntu-2004.dockerfile --build-arg baseimage=nvidia/cuda:11.2.2-devel-ubuntu20.04 --tag eschnett/spack-ubuntu-2004:gpu .
# docker save -o spack-ubuntu-2004-gpu.tar eschnett/spack-ubuntu-2004:gpu
docker push eschnett/spack-ubuntu-2004:gpu
```

## Create a Docker image containing AMReX and other Cactus dependencies
You don't need to do this; these commands are for maintainers.
```sh
cd utils/Docker/make-image-gpu
docker build --tag eschnett/simfactory4-cactus:gpu .
docker push eschnett/simfactory4-cactus:gpu
```
Note: CUDA is a big library and needs 2 GByte of disk space.

## Use the Docker image to build and run Cactus
This is what you need to do.
```sh
./utils/Docker/make/image-cpu/run-container-cpu
```
This opens a shell running in a container using the image above.
