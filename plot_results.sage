import argparse
import itertools
import os
import pdb


def contacts(n,d): # cd for uniform distribution
    if d % 2 == 1:
        return 0
    return 1.0 / (d/2 +1) * binomial(d, d/2) * binomial(2*n - d -1, n - d/2 -1)


def leaves(n, k):
    return (1.0 / n) * binomial(n, k) * binomial(n, k-1)


def root_deg(n, r):
    return (r / (n * 1.0)) * binomial(2*n -1 - r, n - 1)


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


def make_plots(source_base_name, output_base_name, prefix, expectation_function, args):
        source_name = os.path.join('data', prefix + source_base_name)
        with open(source_name, 'r') as f:
            experimental_data = list(map(int, f.readlines()))
        r = args.num_samples / catalan_number(args.n)
        expectation_data = [r * expectation_function(args.n, x) for x in range(len(experimental_data))]
        reg_plot = construct_plot(expectation_data, experimental_data, use_log=False)
        log_plot = construct_plot(expectation_data, experimental_data, use_log=True)

        for (my_plot, plot_prefix) in ((reg_plot, ''), (log_plot, 'log_')):
            plot_name = prefix + plot_prefix + output_base_name
            plot_path = os.path.join('plots', plot_name)
            print('saving {} to: {}'.format(prefix[:-1], plot_path))
            my_plot.save(plot_path)


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


    source_base_name = 'n={}_dist={}_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(args.n, distribution, args.mixing_time, args.sample_interval, args.num_samples)
    plot_base_name = 'n={}_dist={}.png'.format(args.n, distribution)

    make_plots(source_base_name, plot_base_name, 'cd_sums_', contacts, args)
    make_plots(source_base_name, plot_base_name, 'num_leaves_', leaves, args)
    make_plots(source_base_name, plot_base_name, 'root_degree_', root_deg, args)
    
    # plt.ylabel('frequency')
    # plt.xlabel('contact distance')
    # plt.title('simulation with n={}'.format(n))

    # out_file = open('MCMC_thermo_1000n_50000ini_2000int.txt', 'w')
    # for s in outPutSamples:
    #     out_file.write(str(s)+'\n')       
    # out_file.close()
