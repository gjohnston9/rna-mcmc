import argparse
import os
import pdb


def S(n,d): # cd for uniform distribution
	if d % 2 == 1:
		return 0
	return 1 / (d/2 +1) * binomial(d, d/2) * binomial(2*n - d -1, n - d/2 -1)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('n', type=int, help='value of n to use')
	parser.add_argument('mixing_time', type=int)
	parser.add_argument('sample_interval', type=int)
	parser.add_argument('num_samples', type=int)

	move_selection_group = parser.add_mutually_exclusive_group(required=True)
	move_selection_group.add_argument('--by_count', action='store_true', help='when making a move, switch the ith 1 with the jth 0')
	move_selection_group.add_argument('--by_position', action='store_true', help='when making a move, switch the letter at position i with the letter at position j')

	distribution_selection_group = parser.add_mutually_exclusive_group(required=True)
	distribution_selection_group.add_argument('--uniform', action='store_true', help='use a uniform distribution when choosing whether to make a move')
	distribution_selection_group.add_argument('--nntm', action='store_true', help='use the ratio of energies as predicted by the nearest-neighbor thermodynamic model when choosing whether to make a move')


	args = parser.parse_args()

	if args.by_count:
		move_type = 'by_count'
	else:
		move_type = 'by_position'

	if args.uniform:
		distribution = 'uniform'
	else:
		distribution = 'nntm'

	cd_sums_name = 'sums_n={}_moveType={}_dist={}_mixingTime={}_sampleInterval={}_numOfSamples={}.txt'.format(args.n, move_type, distribution, args.mixing_time, args.sample_interval, args.num_samples)
	cd_sums_name = os.path.join('data', cd_sums_name)

	with open(cd_sums_name, 'r') as f:
		cd_sums = list(map(int, f.readlines()))


	r = args.num_samples / catalan_number(args.n)
	expectations = list_plot([r*S(args.n, d) for d in range(4, args.sample_interval)], color='green', size=5)
	cds_plot = list_plot(cd_sums[4:], color='red', size=5)

	plot_name = 'cds_n={}_moveType={}_dist={}.png'.format(args.n, move_type, distribution)
	plot_path = os.path.join('plots', plot_name)
	print('saving figure to: {}'.format(plot_path))
	(cds_plot + expectations).save(plot_path)

	# plt.ylabel('frequency')
	# plt.xlabel('contact distance')
	# plt.title('simulation with n={}'.format(n))

	# out_file = open('MCMC_thermo_1000n_50000ini_2000int.txt', 'w')
	# for s in outPutSamples:
	#     out_file.write(str(s)+'\n')		
	# out_file.close()
