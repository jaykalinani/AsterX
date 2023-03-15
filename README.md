<img align="top" src="Docs/figures/asterx.png" width="140">

**AsterX** is a GPU-accelerated GRMHD code for dynamical spacetimes, written in C++. It is built upon the [CarpetX](https://github.com/eschnett/CarpetX) driver, which is intended for the [Einstein Toolkit](https://einsteintoolkit.org/). **CarpetX** is based on [AMReX](https://amrex-codes.github.io), a software framework for block-structured AMR (adaptive mesh refinement).

Full documentation will soon be available at [asterx.readthedocs.io](https://asterx.readthedocs.io/en/latest/#).

* [![GitHub CI](https://github.com/jaykalinani/AsterX/workflows/CI/badge.svg)](https://github.com/jaykalinani/AsterX/actions)  [![Documentation Status](https://readthedocs.org/projects/asterx/badge/?version=latest)](https://asterx.readthedocs.io/en/latest/?badge=latest) [![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://github.com/jaykalinani/AsterX/blob/main/LICENSE.md)

## Overview 

* Heavily derived from the GRMHD code [Spritz](https://zenodo.org/record/4350072).
* Solves the GRMHD equations in 3D Cartesian coordinates and on dynamical spacetimes using high-resolution shock capturing (HRSC) schemes. 
* Based on the flux-conservative Valencia formulation.
* Directly evolves the staggered vector potential.

## Available modules 

* `AsterX` - the core GRMHD module
* `AsterSeeds` - initial data module
* `Con2PrimFactory` - module providing different conservative-to-primitive variable recovery routines
* `EOSX` - equation of state driver
* `TOVSolver` - a modified version of the publicly available TOVSolver thorn used within the Einstein Toolkit

## Getting started ##

Instructions for downloading and building the Einstein Toolkit including
CarpetX can be found [here](https://github.com/eschnett/CarpetX).

Details for building and running AsterX along with CarpetX will be added to [asterx.readthedocs.io](https://asterx.readthedocs.io/en/latest/#) soon..
