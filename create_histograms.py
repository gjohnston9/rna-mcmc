import matplotlib.pyplot as plt
import numpy as np

from argparse import Namespace
import os

def create_histogram(characteristic_name, uniform, nntm, save_name):
	fig, ax = plt.subplots()
	ax.hist(uniform, normed=True, color='green', alpha=0.5, label='uniform distribution')
	ax.hist(nntm, normed=True, color='red', alpha=0.5, label='thermodynamic distribution')

	ax.set(
		title='Comparison between {} under\nuniform and thermodynamic distributions'.format(characteristic_name),
		xlabel=characteristic_name,
		ylabel='frequency (normalized)')
	ax.legend()
	plt.savefig(save_name)

if __name__ == '__main__':
	n = 1000
	unif = Namespace(mixingTime=0, sampleInterval=1000, numSamples=10000)
	nntm = Namespace(mixingTime=0, sampleInterval=10000, numSamples=1000)
	unif_data_name = os.path.join('data', 'by_sample',
		'avg_branching_n={}_dist=uniform_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(
			n, unif.mixingTime, unif.sampleInterval, unif.numSamples))
	nntm_data_name = os.path.join('data', 'by_sample',
		'avg_branching_n={}_dist=nntm_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(
			n, nntm.mixingTime, nntm.sampleInterval, nntm.numSamples))
	with open(unif_data_name, 'r') as f:
		unif_data = list(map(lambda line: float(line.strip()), f))
	with open(nntm_data_name, 'r') as f:
		nntm_data = list(map(lambda line: float(line.strip()), f))
	### caution: long filename incoming
	save_name = 'avg_branching_n={}_unif_mixing={}_gap={}_samples={}_nntm_mixing={}_gap={}_samples={}.png'.format(
		n,
		unif.mixingTime, unif.sampleInterval, unif.numSamples,
		nntm.mixingTime, nntm.sampleInterval, nntm.numSamples)
	save_name_path = os.path.join('plots', save_name)
	create_histogram('average branching', unif_data, nntm_data, save_name_path)
