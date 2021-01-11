The repository contains the source code for the numerical experiments considered
in [Nonlinear Conjugate Gradient Methods for PDE Constrained Shape Optimization Based on Steklov-Poincaré-Type Metrics](https://arxiv.org/abs/2007.12891) by Sebastian Blauth.

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

Finally, there are two post processing functionalities available in the "visualization"
folder: `optimization_history.py` creates plots showing the evolution of the cost functional and shape gradient norm over the optimization and saves them as .pdf files, and `performance_analysis.py` computes how many iterations the methods need to reach a certain tolerance, and generates the LaTeX tables used in the manuscript.