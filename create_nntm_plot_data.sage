import argparse
import os

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


def write_plot_data(base_name, stat_prefix, expectation_function, n, num_samples):
    source_name = os.path.join('data', 'by_frequency', stat_prefix + base_name)
    with open(source_name, 'r') as f:
        experimental_data = list(map(int, f.readlines()))
    r = num_samples / catalan_number(n)
    expectation_data = [r * expectation_function(n, x) for x in range(1, len(experimental_data)+1)]
    assert len(expectation_data) == len(experimental_data)

    expected_sum = sum(num.n() for num in expectation_data)
    experimental_sum = sum(experimental_data)
    print('stat: {}\nexpected sum: {}\nexperimental sum: {}'.format(stat_prefix[:-1], expected_sum, experimental_sum))

    output_path = os.path.join('data', 'processed_plot_data', 'plotData_' + stat_prefix + base_name)
    print('saving {} to: {}'.format(stat_prefix[:-1], output_path))
    sep = '\t'
    with open(output_path, 'w') as f:
        f.write('experimental frequency under nntm distribution{}expected frequency under uniform distribution\n'.format(sep))
        for experimental, expected in zip(experimental_data, expectation_data):
            f.write('{}{}{}\n'.format(experimental, sep, expected.n()))


def write_plot_data_experimental(nntm_base_name, uniform_base_name, stat_prefix, sub_dir):
    nntm_source_name = os.path.join('data', sub_dir, stat_prefix + nntm_base_name)
    uniform_source_name = os.path.join('data', sub_dir, stat_prefix + uniform_base_name)
    with open(nntm_source_name, 'r') as f:
        nntm_experimental_data = list(map(float, f.readlines()))
    with open(uniform_source_name, 'r') as f:
        uniform_experimental_data = list(map(float, f.readlines()))
    assert len(nntm_experimental_data) == len(uniform_experimental_data)

    output_path = os.path.join('data', 'processed_plot_data', 'plotData_' + stat_prefix + nntm_base_name)
    print('saving {} to: {}'.format(stat_prefix[:-1], output_path))
    sep = '\t'
    with open(output_path, 'w') as f:
        f.write('experimental frequency under nntm distribution{}experimental frequency under uniform distribution\n'.format(sep))
        for nntm, unif in zip(nntm_experimental_data, uniform_experimental_data):
            f.write('{}{}{}\n'.format(nntm, sep, unif))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('n', type=int, help='value of n to use')

    parser.add_argument('uniform_mixing_time', type=int)
    parser.add_argument('uniform_sample_interval', type=int)
    parser.add_argument('uniform_num_samples', type=int)

    parser.add_argument('nntm_params', nargs='*', action='store', help='mixing time, sample interval, and number of samples for nntm distribution')

    distribution_selection_group = parser.add_mutually_exclusive_group(required=True)
    distribution_selection_group.add_argument('--uniform', action='store_true')
    distribution_selection_group.add_argument('--nntm', action='store_true')

    args = parser.parse_args()

    if args.nntm and len(args.nntm_params) == 0:
        error_string = ('If producing data for nntm distribution, provide both nntm and uniform parameters. ' +
            'For characteristics with no available CDF under the uniform distribution, data from a uniform run ' +
            'with the provided uniform parameters will be used instead.')
        parser.error(error_string)

    if args.uniform and len(args.nntm_params) != 0:
        parser.error('If producing data for uniform distribution, do not provide nntm parameters.')

    if args.nntm:
        args.nntm_mixing_time, args.nntm_sample_interval, args.nntm_num_samples = map(int, args.nntm_params)
        nntm_input_name = 'n={}_dist=nntm_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(
            args.n, args.nntm_mixing_time, args.nntm_sample_interval, args.nntm_num_samples)

    unif_input_name = 'n={}_dist=uniform_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(
        args.n, args.uniform_mixing_time, args.uniform_sample_interval, args.uniform_num_samples)

    min_heights = {h : min_height(args.n, h) for h in range(1, args.n+4)}
    def height(n, h, min_heights=min_heights):
        return min_heights[h] - min_heights[h+1]

    if (args.nntm):
        write_plot_data(nntm_input_name, 'height_', height, args.n, args.nntm_num_samples)
        write_plot_data(nntm_input_name, 'num_leaves_', leaves, args.n, args.nntm_num_samples)
        write_plot_data(nntm_input_name, 'root_degree_', root_deg, args.n, args.nntm_num_samples)
        write_plot_data(nntm_input_name, 'cd_sums_', contacts, args.n, args.nntm_num_samples)
        write_plot_data_experimental(nntm_input_name, unif_input_name, 'ladder_distance_', 'by_frequency')
        write_plot_data_experimental(nntm_input_name, unif_input_name, 'cd_averages_', 'by_sample')
    else:
        write_plot_data(unif_input_name, 'height_', height, args.n, args.uniform_num_samples)
        write_plot_data(unif_input_name, 'num_leaves_', leaves, args.n, args.uniform_num_samples)
        write_plot_data(unif_input_name, 'root_degree_', root_deg, args.n, args.uniform_num_samples)
        write_plot_data(unif_input_name, 'cd_sums_', contacts, args.n, args.uniform_num_samples)
