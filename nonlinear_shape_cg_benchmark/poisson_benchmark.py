# Copyright (C) 2020-2021 Sebastian Blauth

import configparser

from fenics import *

from BenchmarkUtils import generate_cg_benchmark_configs, save_results
from cashocs import ShapeOptimizationProblem



set_log_level(LogLevel.CRITICAL)

config = configparser.ConfigParser()
config.read('./poisson_config.ini')
config_list = generate_cg_benchmark_configs(config)

meshlevel = 50
degree = 1
dim = 2
mesh = UnitDiscMesh.create(MPI.comm_world, meshlevel, degree, dim)
coordinates = mesh.coordinates().copy()
dx = Measure('dx', mesh)
ds = Measure('ds', mesh)

boundary = CompiledSubDomain('on_boundary')
boundaries = MeshFunction('size_t', mesh, dim=1)
boundary.mark(boundaries, 1)

V = FunctionSpace(mesh, 'CG', 1)

bcs = DirichletBC(V, Constant(0), boundaries, 1)

x = SpatialCoordinate(mesh)
f = 2.5*pow(x[0] + 0.4 - pow(x[1], 2), 2) + pow(x[0], 2) + pow(x[1], 2) - 1

u = Function(V)
p = Function(V)

e = inner(grad(u), grad(p))*dx - f*p*dx

J = u*dx

for config in config_list:
	mesh.coordinates()[:, :] = coordinates
	optimization_problem = ShapeOptimizationProblem(e, bcs, J, u, p, boundaries, config)
	optimization_problem.solve()

	save_results(config)
