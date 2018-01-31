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


def min_height(n, h):
    """
    returns number of trees of size n of size >= h
    """
    total = 0
    for k in range(1, ((n+1)/h) + 2):
        total += binomial(2*n, n+1-k*h) - 2*binomial(2*n, n-k*h) + binomial(2*n, n-1-k*h)
    return total



def construct_plot(expectation_data, experimental_data, size, xlabel, use_log, truncate_left=False, truncate_right=False):
    if truncate_right:
        for i, item in enumerate(experimental_data[::-1]):
            if item != 0:
                right_lim = len(experimental_data) - i
                break
        else:
            raise ValueError('no non-zero elements found')

        right_lim += 5
        expectation_data = expectation_data[:right_lim]
        experimental_data = experimental_data[:right_lim]

    if truncate_left:
        for i, item in enumerate(experimental_data):
            if item != 0:
                left_lim = i
                break
        else:
            raise ValueError('no non-zero elements found')

        left_lim -= 5
        expectation_data = expectation_data[left_lim:]
        experimental_data = experimental_data[left_lim:]

    if use_log:
        expectation_data = [log(d, base=10) for d in expectation_data]
        experimental_data = [log(d, base=10) for d in experimental_data]

    # to start plotting at x=1
    experimental_datapoints = zip(range(1, len(experimental_data)+1), experimental_data)
    expectation_datapoints = zip(range(1, len(expectation_data)+1), expectation_data)

    plot_expectation = list_plot(expectation_datapoints, color='green', size=size, legend_label='expected values')
    plot_experimental = list_plot(experimental_datapoints, color='red', size=size, legend_label='experimental values')

    max_val = max(itertools.chain(expectation_data, experimental_data))
    min_val = 0

    my_plot = plot_expectation + plot_experimental
    my_plot.set_axes_range(ymin=min_val, ymax=max_val)
    my_plot.axes_labels([xlabel, 'frequency'])
    my_plot.set_legend_options(markerscale=8.0/size)

    return my_plot


def make_plots(source_base_name, output_base_name, prefix, expectation_function, args, size, xlabel, truncate_left=False, truncate_right=False):
        source_name = os.path.join('data', 'by_frequency', prefix + source_base_name)
        with open(source_name, 'r') as f:
            experimental_data = list(map(int, f.readlines()))
        r = args.num_samples / catalan_number(args.n)
        expectation_data = [r * expectation_function(args.n, x) for x in range(1, len(experimental_data))]
        reg_plot = construct_plot(expectation_data, experimental_data, size, xlabel, use_log=False, truncate_left=truncate_left, truncate_right=truncate_right)
        log_plot = construct_plot(expectation_data, experimental_data, size, xlabel, use_log=True, truncate_left=truncate_left, truncate_right=truncate_right)

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

    min_heights = {h : min_height(args.n, h) for h in range(1, args.n+4)}
    def height(n, h, min_heights=min_heights):
        return min_heights[h] - min_heights[h+1]

    make_plots(source_base_name, plot_base_name, 'height_', height, args, 5, 'height', truncate_right=True)
    make_plots(source_base_name, plot_base_name, 'cd_sums_', contacts, args, 5, 'contact distance')
    make_plots(source_base_name, plot_base_name, 'num_leaves_', leaves, args, 5, 'number of leaves', truncate_right=True, truncate_left=True)
    make_plots(source_base_name, plot_base_name, 'root_degree_', root_deg, args, 20, 'root degree', truncate_right=True)

    # plt.ylabel('frequency')
    # plt.xlabel('contact distance')
    # plt.title('simulation with n={}'.format(n))

    # out_file = open('MCMC_thermo_1000n_50000ini_2000int.txt', 'w')
    # for s in outPutSamples:
    #     out_file.write(str(s)+'\n')       
    # out_file.close()
