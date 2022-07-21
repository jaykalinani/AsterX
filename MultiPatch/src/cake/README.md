# Description

The "Cake" patch system consists of 6 Thornburg-like spherical patches that perfectly match with a central cartesian cube patch, thus having a total of 7 patches.

# Patch naming/indexing convention

Each patch in the cake has a name and index. These are defined [here](https://bitbucket.org/eschnett/cactusamrex/src/d89fc96eaf7dbdf831fc76c0fd3a9df8633ee445/MultiPatch/src/cake/cake.hxx#lines-80:96).

# References

The equations that define this patch system can be found in Sec. IV of https://arxiv.org/pdf/gr-qc/0512001.pdf. Other useful references are

1. https://arxiv.org/abs/gr-qc/0404059
2. https://arxiv.org/abs/gr-qc/0507004
3. https://arxiv.org/abs/gr-qc/0512001
4. https://arxiv.org/abs/gr-qc/0602104
5. https://arxiv.org/abs/0712.0353
6. https://arxiv.org/abs/1212.1191
7. https://gwic.ligo.org/assets/docs/theses/reisswig_thesis.pdf
8. https://arxiv.org/abs/0904.0493

# Mathematica Notebooks

Mathematica Notebooks used for visualization and code generation of patches can be found in the [`resources`](resources/) folder.

## Visualizing the Cake patch system

The notebook [`Cake7_Vis.nb`](resources/Cake7_Vis.nb) allows one to interactively visualize the Cake patch system in 3D by drawing constant `c` coordinate surfaces. To produce and interactive plot, simply evaluate the whole notebook.

## Generating C++ code

Jacobians and their derivatives are auto-generated in this patch system. To generate the source files, evaluate the whole notebook file [`Cake7_Gen.nb`](resources/Cake_Gen.nb)

# Pluto.jl notebook

The `Pluto.jl` notebook [`coordinate_tests.jl`](./resources/coordinate_tests.jl) compares data output by `CarpetX` via the [`run_cake_tests.par`](../../par/run_cake_tests.par) parameter file with the expected output of `global2local` and `local2global`. If the tests pass, no output is produced, otherwise an error message will be printed. Note that in order to reproduce the tests, one must change the directory information near the top of the notebook