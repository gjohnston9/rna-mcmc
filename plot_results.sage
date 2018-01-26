import argparse
import itertools
import os
import pdb


def S(n,d): # cd for uniform distribution
    if d % 2 == 1:
        return 0
    return 1 / (d/2 +1) * binomial(d, d/2) * binomial(2*n - d -1, n - d/2 -1)


def construct_plot(data1, data2, use_log):
    if use_log:
        data1 = [log(d, base=10) for d in data1]
        data2 = [log(d, base=10) for d in data2]

    plot1 = list_plot(data1, color='green', size=5)
    plot2 = list_plot(data2, color='red', size=5)

    max_val = max(itertools.chain(data1, data2))
    min_val = 0

    my_plot = plot1 + plot2
    my_plot.set_axes_range(ymin=min_val, ymax=max_val)

    return my_plot


def solo_plot(data, use_log):
    if use_log:
        data = [log(d, base=10) for d in data]

    my_plot = list_plot(data, color='red', size=5)

    max_val = max(data)
    min_val = 0

    my_plot.set_axes_range(ymin=min_val, ymax=max_val)

    return my_plot


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('n', type=int, help='value of n to use')
    parser.add_argument('mixing_time', type=int)
    parser.add_argument('sample_interval', type=int)
    parser.add_argument('num_samples', type=int)

    distribution_selection_group = parser.add_mutually_exclusive_group(required=True)
    distribution_selection_group.add_argument('--uniform', action='store_true', help='use a uniform distribution when choosing whether to make a move')
    distribution_selection_group.add_argument('--nntm', action='store_true', help='use the ratio of energies as predicted by the nearest-neighbor thermodynamic model when choosing whether to make a move')

    args = parser.parse_args()

    if args.uniform:
        distribution = 'uniform'
    else:
        distribution = 'nntm'
    
    cd_sums_name = 'cd_sums_n={}_dist={}_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(args.n, distribution, args.mixing_time, args.sample_interval, args.num_samples)
    cd_sums_name = os.path.join('data', cd_sums_name)

    with open(cd_sums_name, 'r') as f:
        cd_sums_data = list(map(int, f.readlines()))

    r = args.num_samples / catalan_number(args.n)
    expectation_data = [r*S(args.n, d) for d in range(4, len(cd_sums_data))]

    cd_sums_reg_plot = construct_plot(expectation_data, cd_sums_data[4:], use_log=False)
    cd_sums_log_plot = construct_plot(expectation_data, cd_sums_data[4:], use_log=True)





    num_leaves_name = 'num_leaves_n={}_dist={}_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(args.n, distribution, args.mixing_time, args.sample_interval, args.num_samples)
    num_leaves_name = os.path.join('data', num_leaves_name)

    with open(num_leaves_name, 'r') as f:
        num_leaves_data = list(map(int, f.readlines()))

    num_leaves_reg_plot = solo_plot(num_leaves_data, use_log=False)
    num_leaves_log_plot = solo_plot(num_leaves_data, use_log=True)





    base_plot_name = 'n={}_dist={}.png'.format(args.n, distribution)
    
    for (my_plot, prefix, name) in (
        (cd_sums_reg_plot, 'cd_sums_', 'cd_sums_reg_plot'),
        (cd_sums_log_plot, 'cd_sums_log_', 'cd_sums_log_plot'),
        (num_leaves_reg_plot, 'num_leaves_', 'num_leaves_reg_plot'),
        (num_leaves_log_plot, 'num_leaves_log_', 'num_leaves_log_plot')):

        plot_name = prefix + base_plot_name
        plot_path = os.path.join('plots', plot_name)
        print('saving {} to: {}'.format(name, plot_path))
        my_plot.save(plot_path)

    # plt.ylabel('frequency')
    # plt.xlabel('contact distance')
    # plt.title('simulation with n={}'.format(n))

    # out_file = open('MCMC_thermo_1000n_50000ini_2000int.txt', 'w')
    # for s in outPutSamples:
    #     out_file.write(str(s)+'\n')       
    # out_file.close()
