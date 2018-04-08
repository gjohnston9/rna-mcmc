import matplotlib.pyplot as plt
import numpy as np

from argparse import Namespace
import math
import os


def create_histogram(characteristic_name, characteristic_prefix, n, unif_params, nntm_params):
	print('creating histogram for {}'.format(characteristic_name))

	unif_mixingTime, unif_sampleInterval, unif_numSamples = unif_params
	nntm_mixingTime, nntm_sampleInterval, nntm_numSamples = nntm_params

	unif_data_name = '{}n={}_dist=uniform_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(
			characteristic_prefix, n, unif_mixingTime, unif_sampleInterval, unif_numSamples)

	nntm_data_name = '{}n={}_dist=nntm_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(
			characteristic_prefix, n, nntm_mixingTime, nntm_sampleInterval, nntm_numSamples)

	unif_data_name = os.path.join('data', 'by_sample', unif_data_name)
	nntm_data_name = os.path.join('data', 'by_sample', nntm_data_name)

	with open(unif_data_name, 'r') as f:
		uniform = list(map(lambda line: float(line.strip()), f))

	with open(nntm_data_name, 'r') as f:
		nntm = list(map(lambda line: float(line.strip()), f))

	save_name = '{}log_n={}_unif_mixing={}_gap={}_samples={}_nntm_mixing={}_gap={}_samples={}.png'.format(
		characteristic_prefix,
		n,
		unif_mixingTime, unif_sampleInterval, unif_numSamples,
		nntm_mixingTime, nntm_sampleInterval, nntm_numSamples)
	save_name = os.path.join('plots', save_name)


	uniform = [math.log(x) for x in uniform]
	nntm = [math.log(x) for x in nntm]

	min_x = min(min(uniform), min(nntm))
	max_x = max(max(uniform), max(nntm))
	print(max_x)
	padding = (max_x - min_x) / 8.0 # adding an extra 1/8 of whitespace before and after data
	min_x -= padding
	max_x += padding
	min_x = max(0, min_x)
	bins = np.linspace(min_x, max_x, 76)

	fig, ax = plt.subplots()
	ax.hist(nntm, bins=bins, color='red', alpha=0.5, label='thermodynamic distribution')
	ax.hist(uniform, bins=bins, color='green', alpha=0.5, label='uniform distribution')

	ax.set(
		title='Comparison between {} under\nuniform and thermodynamic distributions'.format(characteristic_name),
		xlabel=characteristic_name,
		ylabel='log of frequency')
	ax.legend(loc='upper center', bbox_to_anchor=(0., -0.22, 1., .102), fancybox=False, shadow=False)
	plt.gcf().subplots_adjust(bottom=0.22)
	print('saving {} to {}'.format(characteristic_name, save_name))
	plt.savefig(save_name)


if __name__ == '__main__':
	n = 1000

	unif = Namespace(mixingTime=100000, sampleInterval=1000, numSamples=10000)
	nntm = Namespace(mixingTime=10000000, sampleInterval=10000, numSamples=10000)

	create_histogram(
		'average branching',
		'avg_branching_',
		n,
		(unif.mixingTime, unif.sampleInterval, unif.numSamples),
		(nntm.mixingTime, nntm.sampleInterval, nntm.numSamples))

	create_histogram(
		'contact distance average',
		'cd_averages_',
		n,
		(unif.mixingTime, unif.sampleInterval, unif.numSamples),
		(nntm.mixingTime, nntm.sampleInterval, nntm.numSamples))

