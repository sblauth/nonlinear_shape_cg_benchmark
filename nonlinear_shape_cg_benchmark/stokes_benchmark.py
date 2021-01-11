# Copyright (C) 2020-2021 Sebastian Blauth

import configparser

from fenics import *

from BenchmarkUtils import generate_cg_benchmark_configs, save_results
from cashocs import ShapeOptimizationProblem, import_mesh



set_log_level(LogLevel.CRITICAL)

config = configparser.ConfigParser()
config.read('./stokes_config.ini')
config_list = generate_cg_benchmark_configs(config)

mesh, subdomains, boundaries, dx, ds, dS = import_mesh('./meshes/stokes/mesh.xdmf')
coordinates = mesh.coordinates().copy()

v_elem = VectorElement('CG', mesh.ufl_cell(), 2)
p_elem = FiniteElement('CG', mesh.ufl_cell(), 1)
space = FunctionSpace(mesh, MixedElement([v_elem, p_elem]))

v_in = Expression(('-1.0/4.0*(x[1] - 2.0)*(x[1] + 2.0)', '0.0'), degree=2)
bc_in = DirichletBC(space.sub(0), v_in, boundaries, 1)
bc_wall = DirichletBC(space.sub(0), Constant((0, 0)), boundaries, 2)
bc_gamma = DirichletBC(space.sub(0), Constant((0, 0)), boundaries, 4)
bcs = [bc_in, bc_wall, bc_gamma]


up = Function(space)
u, p = split(up)
vq = Function(space)
v, q = split(vq)

e = inner(grad(u), grad(v))*dx - p*div(v)*dx - q*div(u)*dx

J = inner(grad(u), grad(u))*dx


for config in config_list:
	mesh.coordinates()[:, :] = coordinates
	optimization_problem = ShapeOptimizationProblem(e, bcs, J, up, vq, boundaries, config)
	optimization_problem.solve()

	save_results(config)
