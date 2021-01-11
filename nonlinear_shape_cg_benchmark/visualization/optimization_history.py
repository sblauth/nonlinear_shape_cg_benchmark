# Copyright (C) 2020-2021 Sebastian Blauth

import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np



plt.rcParams.update({'font.size': 14})
lw = 1
ms = 3
alpha = 0.25

problem_list = ['poisson', 'eit', 'stokes', 'pipe']

plt.style.use('tableau-colorblind10')

close_figs = False
save_figs = True

for problem in problem_list:

	with open('../results/' + problem + '/history_gd.json') as file:
		history_gd = json.load(file)

	with open('../results/' + problem + '/history_lbfgs_1.json') as file:
		history_lbfgs_1 = json.load(file)
	with open('../results/' + problem + '/history_lbfgs_3.json') as file:
		history_lbfgs_3 = json.load(file)
	with open('../results/' + problem + '/history_lbfgs_5.json') as file:
		history_lbfgs_5 = json.load(file)

	with open('../results/' + problem + '/history_cg_FR.json') as file:
		history_cg_FR = json.load(file)
	with open('../results/' + problem + '/history_cg_PR.json') as file:
		history_cg_PR = json.load(file)
	with open('../results/' + problem + '/history_cg_HS.json') as file:
		history_cg_HS = json.load(file)
	with open('../results/' + problem + '/history_cg_DY.json') as file:
		history_cg_DY = json.load(file)
	with open('../results/' + problem + '/history_cg_HZ.json') as file:
		history_cg_HZ = json.load(file)

	histories = [history_cg_FR, history_cg_PR, history_cg_HS, history_cg_DY, history_cg_HZ]

	# determine the best cg method
	alphas = [alpha]*5
	iterations = [history['iterations'] for history in histories]
	idx = np.argsort(iterations)
	alphas[idx[0]] = 1

	fig, ax = plt.subplots(figsize=(6, 3))
	if problem == 'poisson':
		ax.plot(range(history_gd['iterations']+1), np.array(history_gd['cost_function_value'])+0.1, label='GD', lw=lw, marker='s', markersize=ms, ls=':', alpha=1)

		# ax.plot(range(history_lbfgs_1['iterations']+1), np.array(history_lbfgs_1['cost_function_value'])+0.1, label='LBFGS 1', lw=lw, marker='^', markersize=ms, ls=':', alpha=alpha)
		# ax.plot(range(history_lbfgs_3['iterations']+1), np.array(history_lbfgs_3['cost_function_value'])+0.1, label='LBFGS 3', lw=lw, marker='^', markersize=ms, ls=':', alpha=alpha)
		ax.plot(range(history_lbfgs_5['iterations']+1), np.array(history_lbfgs_5['cost_function_value'])+0.1, label='LBFGS 5', lw=lw, marker='^', markersize=ms, ls=':', alpha=1)

		ax.plot(range(history_cg_FR['iterations']+1), np.array(history_cg_FR['cost_function_value'])+0.1, label='NCG FR', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[0])
		ax.plot(range(history_cg_PR['iterations']+1), np.array(history_cg_PR['cost_function_value'])+0.1, label='NCG PR', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[1])
		ax.plot(range(history_cg_HS['iterations']+1), np.array(history_cg_HS['cost_function_value'])+0.1, label='NCG HS', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[2])
		ax.plot(range(history_cg_DY['iterations']+1), np.array(history_cg_DY['cost_function_value'])+0.1, label='NCG DY', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[3])
		ax.plot(range(history_cg_HZ['iterations']+1), np.array(history_cg_HZ['cost_function_value'])+0.1, label='NCG HZ', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[4])

	elif problem == 'stokes':
		ax.plot(range(101), history_gd['cost_function_value'][:101], label='GD', lw=lw, marker='s', markersize=ms, ls=':', alpha=1)

		# ax.plot(range(101), history_lbfgs_1['cost_function_value'][:101], label='LBFGS 1', lw=lw, marker='^', markersize=ms, ls=':', alpha=alpha)
		# ax.plot(range(101), history_lbfgs_3['cost_function_value'][:101], label='LBFGS 3', lw=lw, marker='^', markersize=ms, ls=':', alpha=alpha)
		ax.plot(range(history_lbfgs_5['iterations']+1), history_lbfgs_5['cost_function_value'], label='LBFGS 5', lw=lw, marker='^', markersize=ms, ls=':', alpha=1)

		ax.plot(range(101), history_cg_FR['cost_function_value'][:101], label='NCG FR', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[0])
		ax.plot(range(101), history_cg_PR['cost_function_value'][:101], label='NCG PR', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[1])
		ax.plot(range(101), history_cg_HS['cost_function_value'][:101], label='NCG HS', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[2])
		ax.plot(range(history_cg_DY['iterations']+1), history_cg_DY['cost_function_value'], label='NCG DY', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[3])
		ax.plot(range(101), history_cg_HZ['cost_function_value'][:101], label='NCG HZ', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[4])

	else:
		ax.plot(range(history_gd['iterations']+1), history_gd['cost_function_value'], label='GD', lw=lw, marker='s', markersize=ms, ls=':', alpha=1)

		# ax.plot(range(history_lbfgs_1['iterations']+1), history_lbfgs_1['cost_function_value'], label='LBFGS 1', lw=lw, marker='^', markersize=ms, ls=':', alpha=alpha)
		# ax.plot(range(history_lbfgs_3['iterations']+1), history_lbfgs_3['cost_function_value'], label='LBFGS 3', lw=lw, marker='^', markersize=ms, ls=':', alpha=alpha)
		ax.plot(range(history_lbfgs_5['iterations']+1), history_lbfgs_5['cost_function_value'], label='LBFGS 5', lw=lw, marker='^', markersize=ms, ls=':', alpha=1)

		ax.plot(range(history_cg_FR['iterations']+1), history_cg_FR['cost_function_value'], label='NCG FR', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[0])
		ax.plot(range(history_cg_PR['iterations']+1), history_cg_PR['cost_function_value'], label='NCG PR', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[1])
		ax.plot(range(history_cg_HS['iterations']+1), history_cg_HS['cost_function_value'], label='NCG HS', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[2])
		ax.plot(range(history_cg_DY['iterations']+1), history_cg_DY['cost_function_value'], label='NCG DY', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[3])
		ax.plot(range(history_cg_HZ['iterations']+1), history_cg_HZ['cost_function_value'], label='NCG HZ', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[4])

	if problem == 'stokes':
		ax.set_xticks(np.arange(0, 101, 20))
	else:
		ax.set_xticks(np.arange(0, history_gd['iterations'] + 1, 10))
	ax.legend(prop={'size': 10})
	ax.set_xlabel('Iterations')
	ax.set_ylabel('Cost Functional Value')
	if problem == 'eit' or problem=='poisson':
		ax.set_yscale('log')
		ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=10))
		ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), numticks=10))
	elif problem == 'pipe':
		ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
	if save_figs:
		fig.savefig('./pdfs/cost_functional_' + problem + '.pdf', dpi=1000, bbox_inches='tight')
	if close_figs:
		plt.close(fig)



	fig, ax = plt.subplots(figsize=(6, 3))
	ax.plot(range(history_gd['iterations']+1), history_gd['gradient_norm'], label='GD', lw=lw, marker='s', markersize=ms, ls=':', alpha=1)

	# ax.plot(range(history_lbfgs_1['iterations']+1), history_lbfgs_1['gradient_norm'], label='LBFGS 1', lw=lw, marker='^', markersize=ms, ls=':', alpha=alpha)
	# ax.plot(range(history_lbfgs_3['iterations']+1), history_lbfgs_3['gradient_norm'], label='LBFGS 3', lw=lw, marker='^', markersize=ms, ls=':', alpha=alpha)
	ax.plot(range(history_lbfgs_5['iterations']+1), history_lbfgs_5['gradient_norm'], label='LBFGS 5', lw=lw, marker='^', markersize=ms, ls=':', alpha=1)

	ax.plot(range(history_cg_FR['iterations']+1), history_cg_FR['gradient_norm'], label='NCG FR', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[0])
	ax.plot(range(history_cg_PR['iterations']+1), history_cg_PR['gradient_norm'], label='NCG PR', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[1])
	ax.plot(range(history_cg_HS['iterations']+1), history_cg_HS['gradient_norm'], label='NCG HS', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[2])
	ax.plot(range(history_cg_DY['iterations']+1), history_cg_DY['gradient_norm'], label='NCG DY', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[3])
	ax.plot(range(history_cg_HZ['iterations']+1), history_cg_HZ['gradient_norm'], label='NCG HZ', lw=lw, marker='o', markersize=ms, ls='-', alpha=alphas[4])

	if problem == 'stokes':
		ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(50))
	else:
		ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))

	
	ax.legend(prop={'size': 10})
	ax.set_yscale('log')
	ax.set_xlabel('Iterations')
	
	if problem == 'eit':
		ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=10))
		ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), numticks=10))
	ax.set_ylabel('Shape Gradient Norm')
	lines = [pow(10, -i) for i in range(1,int(-np.ceil(np.log10(ax.get_ylim()[0])))+1)]
	ax.hlines(lines, 0, history_gd['iterations'], ls=':')
	if save_figs:
		fig.savefig('./pdfs/gradient_norm_' + problem + '.pdf', dpi=1000, bbox_inches='tight')
	if close_figs:
		plt.close(fig)
