# Copyright (C) 2020-2021 Sebastian Blauth

import copy
import os
import shutil



def generate_cg_benchmark_configs(config):
	"""

	Parameters
	----------
	config : configparser.ConfigParser

	Returns
	-------

	"""
	
	config_list = []

	### Gradient Descent
	temp_config = copy.deepcopy(config)
	temp_config.set('OptimizationRoutine', 'algorithm', 'gd')
	config_list.append(temp_config)

	### BFGS
	temp_config = copy.deepcopy(config)
	temp_config.set('OptimizationRoutine', 'algorithm', 'bfgs')
	temp_config.set('AlgoLBFGS', 'bfgs_memory_size', '1')
	config_list.append(temp_config)

	temp_config = copy.deepcopy(config)
	temp_config.set('OptimizationRoutine', 'algorithm', 'bfgs')
	temp_config.set('AlgoLBFGS', 'bfgs_memory_size', '3')
	config_list.append(temp_config)

	temp_config = copy.deepcopy(config)
	temp_config.set('OptimizationRoutine', 'algorithm', 'bfgs')
	temp_config.set('AlgoLBFGS', 'bfgs_memory_size', '5')
	config_list.append(temp_config)

	### CG
	temp_config = copy.deepcopy(config)
	temp_config.set('OptimizationRoutine', 'algorithm', 'cg')
	temp_config.set('AlgoCG', 'cg_method', 'FR')
	config_list.append(temp_config)

	temp_config = copy.deepcopy(config)
	temp_config.set('OptimizationRoutine', 'algorithm', 'cg')
	temp_config.set('AlgoCG', 'cg_method', 'PR')
	config_list.append(temp_config)

	temp_config = copy.deepcopy(config)
	temp_config.set('OptimizationRoutine', 'algorithm', 'cg')
	temp_config.set('AlgoCG', 'cg_method', 'HS')
	config_list.append(temp_config)

	temp_config = copy.deepcopy(config)
	temp_config.set('OptimizationRoutine', 'algorithm', 'cg')
	temp_config.set('AlgoCG', 'cg_method', 'DY')
	config_list.append(temp_config)

	temp_config = copy.deepcopy(config)
	temp_config.set('OptimizationRoutine', 'algorithm', 'cg')
	temp_config.set('AlgoCG', 'cg_method', 'HZ')
	config_list.append(temp_config)

	return config_list



def save_results(config):
	
	result_dir = config.get('Output', 'result_dir')
	
	if config.get('OptimizationRoutine', 'algorithm') in ['gd', 'gradient_descent']:
		shutil.rmtree(result_dir + '/pvd_gd', ignore_errors=True)
		try:
			os.remove(result_dir + '/history_gd.json')
		except:
			pass

		try:
			os.rename(result_dir + '/pvd', result_dir + '/pvd_gd')
		except:
			pass
		os.rename(result_dir + '/history.json', result_dir + '/history_gd.json')


	elif config.get('OptimizationRoutine', 'algorithm') in ['bfgs', 'lbfgs']:
		shutil.rmtree(result_dir + '/pvd_lbfgs_' + config.get('AlgoLBFGS', 'bfgs_memory_size'), ignore_errors=True)
		try:
			os.remove(result_dir + '/history_lbfgs_' + config.get('AlgoLBFGS', 'bfgs_memory_size') + '.json')
		except:
			pass

		try:
			os.rename(result_dir + '/pvd', result_dir + '/pvd_lbfgs_' + config.get('AlgoLBFGS', 'bfgs_memory_size'))
		except:
			pass
		os.rename(result_dir + '/history.json', result_dir + '/history_lbfgs_' + config.get('AlgoLBFGS', 'bfgs_memory_size') + '.json')


	elif config.get('OptimizationRoutine', 'algorithm') in ['conjugate_gradient', 'cg']:
		shutil.rmtree(result_dir + '/pvd_cg_' + config.get('AlgoCG', 'cg_method', fallback='FR'), ignore_errors=True)
		try:
			os.remove(result_dir + '/history_cg_' + config.get('AlgoCG', 'cg_method', fallback='FR') + '.json')
		except:
			pass

		try:
			os.rename(result_dir + '/pvd', result_dir + '/pvd_cg_' + config.get('AlgoCG', 'cg_method', fallback='FR'))
		except:
			pass
		os.rename(result_dir + '/history.json', result_dir + '/history_cg_' + config.get('AlgoCG', 'cg_method', fallback='FR') + '.json')
