# [CactusAMReX](https://bitbucket.org/eschnett/cactusamrex)

![CarpetX logo](figures/carpetx.png)

**CarpetX** is a [Cactus](https://cactuscode.org/) driver based on [AMReX](https://amrex-codes.github.io), a software framework for block-structured AMR (adaptive mesh refinement). CarpetX is intended for the [Einstein Toolkit](https://einsteintoolkit.org/).

* [Bitbucket](https://bitbucket.org/eschnett/cactusamrex): Source code repository
* [Azure   Pipelines](https://dev.azure.com/schnetter/CactusAMReX/_build): Build Status [![Build Status](https://dev.azure.com/schnetter/CactusAMReX/_apis/build/status/CactusAMReX-CI?branchName=master)](https://dev.azure.com/schnetter/CactusAMReX/_build/latest?definitionId=6&branchName=master)

## Overview

CarpetX is almost ready for production. (The only missing feature is
checkpointing/recovery.) You are welcome to give it a try, to look at
what changes your code might need to benefit from CarpetX's new
features, and to give us feedback.

The recorded talk ["Using CarpetX: A Guide for Early
Adopters"](http://einsteintoolkit.org/seminars/2021_03_18/index.html).
This presentation provides an overview of the current capabilities of
CarpetX and showcases how to write Cactus code using it.

## Getting started

Here are instructions for downloading the Einstein Toolkit including
CarpetX, building, and running an example.

### Download and Setup

Download the Einstein Toolkit, including CarpetX. This will create a
new directory `Cactus` that will contain the code:
```sh
curl -kLO https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
perl GetComponents --parallel https://bitbucket.org/eschnett/cactusamrex/raw/master/manifest/einsteintoolkit-carpetx.th
cd Cactus
```

We're using a Docker image to provide dependencies (including AMReX)
to simplify installation. However, the Einstein Toolkit source code
does not reside in that image; it resides in the `Cactus` directory
which you just created. This leads to the following workflow:
- To download, look at, edit, git add/commit/pull/push etc. the code,
  you use the regular tools you already have installed on your system.
  Docker is not involved in this in any way.
- To build and run, you create an emphemeral (stateless) Docker
  container running Bash, with our Docker image mounted. For
  convenience, there is a shell script for this:
```sh
./repos/cactusamrex/utils/Docker/cpu/run-container
```
The first time you run this script, Docker will download the Docker
image with the dependencies, which might take a few minutes.

Inside this container, the hostname is `carpetx-docker-cpu`, so that
you know you're inside the container. (You exit the container by
exiting the shell, e.g. with the `exit` command.)

Note: While the container (running the Bash shell) is thrown away
after each use, the changes to the file system you make persist, since
your home directory is mounted and thus available in the container.

As usual, the first time you use the Einstein Toolkit, you have to
configure Simfactory for your local machine:
```sh
./simfactory/bin/sim setup
```

### Build

Let's build the toolkit with CarpetX.

The default option list doesn't point to AMReX nor some other
dependencies. We thus specify our own. (I assume this could be fixed
in the thorns handling these external dependencies.)

The default thorn list would build the regular Einstein Toolkit
without CarpetX; we need to specify a particular thorn list which
includes all CarpetX thorns. This thorn list also disables those
thorns that do not yet work with CarpetX.

```sh
./simfactory/bin/sim build --optionlist=repos/cactusamrex/utils/Docker/cpu/carpetx.cfg --thornlist=repos/cactusamrex/utils/Docker/cpu/carpetx.th
```

Of course, once you created your configuration with the command above,
to re-build after changing some code, the command is simply
```sh
./simfactory/bin/sim build
```
as usual.

### Run

Let's run some examples!

Starting slowly, here is a scalar wave:
```sh
./simfactory/bin/sim submit planewave --parfile=arrangements/CarpetX/WaveToyCPU/par/planewave.par --procs=8 --num-threads=8
```
You want to adapt the number of cores (`--procs`) and threads
(`--num-threads`) to your system.

It seems that the default Simfactory setup that was created above
buffers the output, so that there might be long periords of time where
the simulation appears to hang. Use `top` to see whether it is still
running, and check the output directory to see whether it is still
producing output.

We can also run a binary black hole merger:
```sh
./simfactory/bin/sim submit planewave --parfile=arrangements/CarpetX/Z4c/par/qc0.rpar --procs=40 --num-threads=10
```

This setup requires more memory and time. I'm running it with about
200 GByte of memory on 40 cores for 24 hours.
