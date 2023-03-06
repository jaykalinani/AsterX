**AsterX** is a GPU-accelerated GRMHD code for dynamical spacetimes, which is built upon the [CarpetX](https://github.com/eschnett/CarpetX) driver.
**CarpetX** is based on [AMReX](https://amrex-codes.github.io), a software framework for block-structured AMR (adaptive mesh refinement), which is intended for the [Einstein Toolkit](https://einsteintoolkit.org/).

* [![GitHub
  CI](https://github.com/jaykalinani/AsterX/workflows/ci/badge.svg)](https://github.com/jaykalinani/AsterX/actions)

## Overview

Most algorithms in AsterX are heavily derived from the GRMHD code [Spritz](https://zenodo.org/record/4350072). Like Spritz, AsterX solves the GRMHD equations in 3D Cartesian coordinates and on a dynamical spacetime using high-resolution shock capturing (HRSC) schemes. AsterX is based on the flux-conservative Valencia formulation, and directly evolves the vector potential with the generalised Lorenz gauge to guarantee that the divergence-free constraint of the magnetic field is always satisfied with very high accuracy (thus avoiding magnetic monopoles). Moreover, the vector potential is defined in a staggered way (like in Spritz and differently from WhiskyMHD).

The flux-conservative approach involves evolving a set of conservative equations for certain conserved quantities. At each time-step, it is required to recover the primitive variables (like density, velocity, and energy density) from conservative ones using particular conservative-to-primitive (C2P) variable conversion numerical schemes. C2P schemes are closely dependent on the equation of state (EOS) of the fluid, and are generally prone to errors/failures in GRMHD simulations (see Siegel et al. 2018 for various implementations of recovery schemes). In AsterX, we have currently implemented the C2P scheme of Noble et al. 2006 which is based on a Newton-Raphson method and the C2P scheme of Palenzuela et al. 2015 which instead is based on the root-bracketing method. Morever, we also take advantage of the new robust and accurate C2P scheme of Kastaun et al. 2021 (which has been tested rigorously in Kalinani et al. 2022), via its publicly available library called RePrimAnd. However, this library is currently supported only on CPUs, but efforts are on the way to make it GPU-compatible as well. All the three C2P schemes are implemented in an EOS agnostic way, however they use an EOS framework which currently supports only analytical EOSs such as polytropic EOS and ideal gas EOS. We are also currently extending this framework to allow for tabulated cold EOS, hybrid EOS (tabulated cold part plus a thermal component) as well as finite temperature, composition dependent tabulated EOS. 

To compute the fluxes, the primitive variables first need to be reconstructed on the cell-interfaces. To do this, we have implemented the second-order accurate total variation diminishing (TVD) 'Minmod' scheme and third-order accurate Piecewise Parabolic reconstruction Method (PPM). To treat the Riemann problems at the cell-interface, we have included second-order accurate Lax-Friedrichs and third-order accurate HLLE flux solvers. For stepping forward in time, we employ ODE solvers such as fourth-order Runge-Kutta methods. The CarpetX infrastructure allows AsterX to exactly conserve the MHD quantities across AMR boundaries through refluxing. At each time-step, space-time is evolved using the Z4c solver which based on a conformal decomposition of the Z4 formulation (Bernuzzi et al. 2010), while the GRMHD equations are evolved via AsterX.

To assess the accuracy and robustness of the basic GRMHD algorithms implemented, AsterX is currently being extensively tested via a series of benchmarks in Minkowski and curved spacetimes which include relativistic Riemann shocktube problems, magnetised cylindrical and spherical explosions, magnetic rotors, magnetic loop advection, as well as evolution of a magnetised oscillating stable neutron star.


## Getting started

Instructions for downloading and building the Einstein Toolkit including
CarpetX can be found [here]((https://github.com/eschnett/CarpetX)).

Details for building and running AsterX along with CarpetX will be added soon..


  
