# Copyright (C) 2020-2021 Sebastian Blauth

import configparser

from fenics import *

from BenchmarkUtils import generate_cg_benchmark_configs, save_results
from cashocs import ShapeOptimizationProblem, import_mesh



set_log_level(LogLevel.CRITICAL)


config = configparser.ConfigParser()
config.read('./pipe_config.ini')
config_list = generate_cg_benchmark_configs(config)

Re = 4e2

mesh, subdomains, boundaries, dx, ds, dS = import_mesh('./meshes/pipe/mesh.xdmf')
coordinates = mesh.coordinates().copy()

v_elem = VectorElement('CG', mesh.ufl_cell(), 2, dim=2)
p_elem = FiniteElement('CG', mesh.ufl_cell(), 1)
space = FunctionSpace(mesh, MixedElement([v_elem, p_elem]))

v_in = Expression(('-6*(x[1] - 1)*(x[1] + 0)', '0.0'), degree=2)
bc_in = DirichletBC(space.sub(0), v_in, boundaries, 1)
bc_wall = DirichletBC(space.sub(0), Constant((0, 0)), boundaries, 2)
bc_gamma = DirichletBC(space.sub(0), Constant((0, 0)), boundaries, 4)
bcs = [bc_in, bc_wall, bc_gamma]


up = Function(space)
u, p = split(up)
vq = Function(space)
v, q = split(vq)

e = Constant(1/Re)*inner(grad(u), grad(v))*dx + inner(grad(u)*u, v)*dx - p*div(v)*dx - q*div(u)*dx

J = Constant(1/Re)*inner(grad(u), grad(u))*dx

initial_guess = [[interpolate(Constant((0,0)), space.sub(0).collapse()),
				 interpolate(Constant(0), space.sub(1).collapse())]]

for config in config_list:
	# if config.get('OptimizationRoutine', 'algorithm') == 'cg' and config.get('AlgoCG', 'cg_method', fallback='FR') == 'DY':
	# 	config.set('AlgoCG', 'cg_periodic_restart', 'true')
	# 	config.set('AlgoCG', 'cg_periodic_its', '10')

	mesh.coordinates()[:, :] = coordinates
	optimization_problem = ShapeOptimizationProblem(e, bcs, J, up, vq, boundaries, config, initial_guess=initial_guess)
	optimization_problem.solve()

	save_results(config)
