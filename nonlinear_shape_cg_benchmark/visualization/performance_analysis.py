# Copyright (C) 2020-2021 Sebastian Blauth

import json

import numpy as np



problem_list = ['poisson', 'eit', 'stokes', 'pipe']

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


	histories = [history_gd, history_lbfgs_1, history_lbfgs_3, history_lbfgs_5, history_cg_FR, history_cg_PR, history_cg_HS, history_cg_DY, history_cg_HZ]
	
	### uncomment for command line output
	# print('Number of PDE solves for ' + problem + ':\n')
	# print('Gradient Descent: ' + str(history_gd['state_solves']) + ' (state)    ' + str(history_gd['adjoint_solves']) + ' (adjoint)')
	#
	# print('LBFGS 1: ' + str(history_lbfgs_1['state_solves']) + ' (state)    ' + str(history_lbfgs_1['adjoint_solves']) + ' (adjoint)')
	# print('LBFGS 3: ' + str(history_lbfgs_3['state_solves']) + ' (state)    ' + str(history_lbfgs_3['adjoint_solves']) + ' (adjoint)')
	# print('LBFGS 5: ' + str(history_lbfgs_5['state_solves']) + ' (state)    ' + str(history_lbfgs_5['adjoint_solves']) + ' (adjoint)')
	#
	# print('CG FR: ' + str(history_cg_FR['state_solves']) + ' (state)    ' + str(history_cg_FR['adjoint_solves']) + ' (adjoint)')
	# print('CG PR: ' + str(history_cg_PR['state_solves']) + ' (state)    ' + str(history_cg_PR['adjoint_solves']) + ' (adjoint)')
	# print('CG HS: ' + str(history_cg_HS['state_solves']) + ' (state)    ' + str(history_cg_HS['adjoint_solves']) + ' (adjoint)')
	# print('CG DY: ' + str(history_cg_DY['state_solves']) + ' (state)    ' + str(history_cg_DY['adjoint_solves']) + ' (adjoint)')
	# print('CG HZ: ' + str(history_cg_HZ['state_solves']) + ' (state)    ' + str(history_cg_HZ['adjoint_solves']) + ' (adjoint)\n')

	minimum = 1e0
	for history in histories:
		temp_min = np.min(history['gradient_norm'])
		if temp_min <= minimum:
			minimum = temp_min


	tol_list = np.array([1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5])
	length = np.max(np.where(tol_list >= minimum)) + 1
	tol_list = tol_list[:length]
	iteration_counter_list = [np.array(['-']*length, dtype=object) for repet in range(9)]

	for idx, iteration_counter in enumerate(iteration_counter_list):
		for tol_idx, tol in enumerate(tol_list):
			try:
				iteration_counter[tol_idx] = str(np.min(np.where(np.array(histories[idx]['gradient_norm']) < tol)))
			except:
				pass


	with open('./txt/performance_results_' + problem + '.txt', 'w') as f:
		strings = ['', '', '', '', '', '', '', '', '']

		for idx, _ in enumerate(strings):
			for i in range(length):
				strings[idx] += '& ' + iteration_counter_list[idx][i] + ' '

			strings[idx] += '&  & ' + str(histories[idx]['state_solves']) + ' / ' + str(histories[idx]['adjoint_solves']) + ' \\\\'

		header = 'tol '
		for i in range(length):
			header += '& \\num{' + format(tol_list[i], '.0e') + '} '
		header += ' &  & state / adjoint solves \\\\'

		print('\\toprule', file=f)
		print(header, file=f)
		print('\\midrule', file=f)
		print('GD ' + strings[0], file=f)
		print('\\midrule', file=f)
		print('L-BFGS 1 ' + strings[1], file=f)
		print('L-BFGS 3 ' + strings[2], file=f)
		print('L-BFGS 5 ' + strings[3], file=f)
		print('\\midrule', file=f)
		print('CG FR ' + strings[4], file=f)
		print('CG PR ' + strings[5], file=f)
		print('CG HS ' + strings[6], file=f)
		print('CG DY ' + strings[7], file=f)
		print('CG HZ ' + strings[8], file=f)
		print('\\bottomrule', file=f)
