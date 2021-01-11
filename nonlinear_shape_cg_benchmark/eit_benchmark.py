# Copyright (C) 2020-2021 Sebastian Blauth

import configparser

from fenics import *

from BenchmarkUtils import generate_cg_benchmark_configs, save_results
from cashocs import ShapeOptimizationProblem, import_mesh



set_log_level(LogLevel.CRITICAL)

sigma_out = 1e0
sigma_in = 1e1

def generate_references():
	mesh, subdomains, boundaries, dx, ds, dS = import_mesh('./meshes/eit/reference.xdmf')

	cg_elem = FiniteElement('CG', mesh.ufl_cell(), 1)
	r_elem = FiniteElement('R', mesh.ufl_cell(), 0)
	V = FunctionSpace(mesh, MixedElement([cg_elem, r_elem]))

	u, c = TrialFunctions(V)
	v, d = TestFunctions(V)

	a = sigma_out*inner(grad(u), grad(v))*dx(1) + sigma_in*inner(grad(u), grad(v))*dx(2) + u*d*ds + v*c*ds
	L1  = Constant(1)*v*(ds(3) + ds(4)) + Constant(-1)*v*(ds(1) + ds(2))
	L2  = Constant(1)*v*(ds(3) + ds(2)) + Constant(-1)*v*(ds(1) + ds(4))
	L3  = Constant(1)*v*(ds(3) + ds(1)) + Constant(-1)*v*(ds(2) + ds(4))

	reference1 = Function(V)
	reference2 = Function(V)
	reference3 = Function(V)
	solve(a==L1, reference1)
	solve(a==L2, reference2)
	solve(a==L3, reference3)

	ref1, _ = reference1.split(True)
	ref2, _ = reference2.split(True)
	ref3, _ = reference3.split(True)

	return [ref1, ref2, ref3]


config = configparser.ConfigParser()
config.read('./eit_config.ini')

config_list = generate_cg_benchmark_configs(config)

mesh, subdomains, boundaries, dx, ds, dS = import_mesh('./meshes/eit/mesh.xdmf')
coordinates = mesh.coordinates().copy()
cg_elem = FiniteElement('CG', mesh.ufl_cell(), 1)
r_elem = FiniteElement('R', mesh.ufl_cell(), 0)
V = FunctionSpace(mesh, MixedElement([cg_elem, r_elem]))

references = generate_references()

bcs = None

uc1 = Function(V)
u1, c1 = split(uc1)
pd1 = Function(V)
p1, d1 = split(pd1)
e1 = sigma_out*inner(grad(u1), grad(p1))*dx(1) + sigma_in*inner(grad(u1), grad(p1))*dx(2) \
	 + u1*d1*ds + p1*c1*ds - Constant(1)*p1*(ds(3) + ds(4)) - Constant(-1)*p1*(ds(1) + ds(2))

uc2 = Function(V)
u2, c2 = split(uc2)
pd2 = Function(V)
p2, d2 = split(pd2)
e2 = sigma_out*inner(grad(u2), grad(p2))*dx(1) + sigma_in*inner(grad(u2), grad(p2))*dx(2) \
	 + u2*d2*ds + p2*c2*ds - Constant(1)*p2*(ds(3) + ds(2)) - Constant(-1)*p2*(ds(1) + ds(4))

uc3 = Function(V)
u3, c3 = split(uc3)
pd3 = Function(V)
p3, d3 = split(pd3)
e3 = sigma_out*inner(grad(u3), grad(p3))*dx(1) + sigma_in*inner(grad(u3), grad(p3))*dx(2) \
	 + u3*d3*ds + p3*c3*ds - Constant(1)*p3*(ds(3) + ds(1)) - Constant(-1)*p3*(ds(2) + ds(4))

e = [e1, e2, e3]
u = [uc1, uc2, uc3]
p = [pd1, pd2, pd3]

mu1 = Expression('val', degree=0, val=1.0, domain=mesh)
mu2 = Expression('val', degree=0, val=1.0, domain=mesh)
mu3 = Expression('val', degree=0, val=1.0, domain=mesh)

J1 = mu1*Constant(0.5)*pow(u1 - references[0], 2)*ds
J2 = mu2*Constant(0.5)*pow(u2 - references[1], 2)*ds
J3 = mu3*Constant(0.5)*pow(u3 - references[2], 2)*ds


J = J1 + J2 + J3

DG0 = FunctionSpace(mesh, 'DG', 0)
post = Function(DG0)
a_post = TrialFunction(DG0)*TestFunction(DG0)*dx
L_post = Constant(1)*TestFunction(DG0)*dx(1) + Constant(2)*TestFunction(DG0)*dx(2)


for config in config_list:
	mesh.coordinates()[:, :] = coordinates

	mu1.val = 1.0
	mu2.val = 1.0
	mu3.val = 1.0

	optimization_problem = ShapeOptimizationProblem(e, bcs, J, u, p, boundaries, config)
	optimization_problem.state_problem.solve()

	mu1.val = 1/assemble(J1)
	mu2.val = 1/assemble(J2)
	mu3.val = 1/assemble(J3)

	optimization_problem.solve()
	save_results(config)

	solve(a_post==L_post, post)

	if config.get('OptimizationRoutine', 'algorithm') in ['gd', 'gradient_descent']:
		post_file = File('./results/eit/postprocessing/gd.pvd')
	elif config.get('OptimizationRoutine', 'algorithm') in ['bfgs', 'lbfgs']:
		post_file = File('./results/eit/postprocessing/lbfgs_' + config.get('AlgoLBFGS', 'bfgs_memory_size') + '.pvd')
	elif config.get('OptimizationRoutine', 'algorithm') in ['cg', 'conjugate_gradient']:
		post_file = File('./results/eit/postprocessing/cg_' + config.get('AlgoCG', 'cg_method', fallback='FR') + '.pvd')

	post_file << post
