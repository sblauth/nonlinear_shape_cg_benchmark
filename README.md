[![DOI](https://img.shields.io/badge/DOI-10.1137%2F20M1367738-blue)](https://doi.org/10.1137/20M1367738)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4432327.svg)](https://doi.org/10.5281/zenodo.4432327)




The repository contains the source code for the numerical experiments considered
in [Nonlinear Conjugate Gradient Methods for PDE Constrained Shape Optimization Based on Steklov-Poincaré-Type Metrics](https://doi.org/10.1137/20M1367738) by Sebastian Blauth.

To run the code, you have to install [cashocs](https://cashocs.readthedocs.io/)
first, which includes all necessary prerequisites. The results presented in this
repository have been obtained with version 1.2.1 of cashocs (which uses FEniCS 2019.1).

The repository consists of the following four test cases:

- A shape optimization problem with a Poisson equation (named `poisson`) which
is considered in Section 4.2 of the manuscript.

- A shape identification problem in electrical impedance tomography (named `eit`) which is considered in Section 4.3 of the manuscript.

- A shape optimization problem in Stokes flow (named `stokes`) which is considered
in Section 4.4 of the manuscript.

- A shape optimization problem for a pipe with Navier-Stokes flow (named `pipe`)
which is considered in Section 4.5 of the manuscript.

For each case, there exist two files: `case_benchmark.py` file, containing
the actual benchmark, and `case_config.ini`, the config file for the respective
setting. Here, `case` is one of `poisson`, `eit`, `stokes`, or `pipe`. In particular,
the benchmark for the Poisson problem can be run with the command

    python3 poisson_benchmark.py

and the other benchmarks can be run analogously.

Note, that each problem is solved by the gradient descent method (abbreviated `gd`), limited memory BFGS methods with a memory size of 1, 3, and 5 (abbreviated `lbfgs_1`, `lbfgs_3`, and `lbfgs_5`), and the Fletcher-Reeves (`cg_FR`), Polak-Ribiere (`cg_PR`),
Hestenes-Stiefel (`cg_HS`), Dai-Yuan (`cg_DY`), and Hager-Zhang (`cg_HZ`) nonlinear
conjugate gradient methods.

Finally, there are two post processing functionalities available in the "visualization"
folder: `optimization_history.py` creates plots showing the evolution of the cost functional and shape gradient norm over the optimization and saves them as .pdf files, and `performance_analysis.py` computes how many iterations the methods need to reach a certain tolerance, and generates the LaTeX tables used in the manuscript.

This software is citeable under the following DOI: [10.5281/zenodo.4432327](https://doi.org/10.5281/zenodo.4432327).

If you use these nonlinear CG methods for your work, please cite the paper

    Nonlinear Conjugate Gradient Methods for PDE Constrained Shape Optimization Based on Steklov-Poincaré-Type Metrics
    Sebastian Blauth
    SIAM Journal on Optimization, Volume 31, Issue 3
    https://doi.org/10.1137/20M1367738

If you are using BibTeX, you can use the following entry:

    @Article{Blauth2020Nonlinear,
    author   = {Sebastian Blauth},
    journal  = {SIAM J. Optim.},
    title    = {{N}onlinear {C}onjugate {G}radient {M}ethods for {PDE} {C}onstrained {S}hape {O}ptimization {B}ased on {S}teklov-{P}oincaré-{T}ype {M}etrics},
    year     = {2021},
    number   = {3},
    pages    = {1658--1689},
    volume   = {31},
    doi      = {10.1137/20M1367738},
    fjournal = {SIAM Journal on Optimization},
    }