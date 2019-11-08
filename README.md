# BELLA package
This matlab package provides generic solvers for Bregman EnveLope Linesearch Algorithm (BELLA) for solving structured nonsmooth nonconvex problems of the form

min_x f(x)+ g(x),

where f is relatively smooth and g is proper and lower semicontinuous.

We provide a solver for solving the above-mentioned problem:

BELLA: Bregman EnveLope Linesearch Algorithm

# Using the package
## Installation
All the functions in the package are given in the folder /MatlabCode.

In order to use the function we recommend to execute the following command

addpath(genpath('.'))

if you are not working in the root folder of the package or replacing '.' by the location of the folder on your machine.

## Example: nonnegative matrix factorization (NMF)
We recommend to look at the following files to see how to use the package:

Demo/demo_NMF.m: contains an example of BELLA for NMF with synthetic data.


# References
[1] M. Ahookhosh, A. Themelis, and P. Patrinos, Bregman forward-backward splitting for nonconvex composite optimization:
superlinear convergence to nonisolated critical points, (2019) https://arxiv.org/pdf/1905.11904.pdf
