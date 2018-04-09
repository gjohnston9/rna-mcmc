import argparse
import itertools
import os
import pdb

def construct_plot(expectation_data, experimental_data, size, xlabel, use_log, dist, truncate_left=False, truncate_right=False):
    assert dist in ('nntm, uniform')
    expectation_data = list(enumerate(expectation_data, 1))
    experimental_data = list(enumerate(experimental_data, 1))

    if truncate_right:
        for i, ((x1, y1), (x2, y2)) in enumerate(zip(experimental_data[::-1], expectation_data[::-1])):
            if y1 >= 1 or y2 >= 1:
                right_lim = len(experimental_data) - i
                break
        else:
            raise ValueError('no non-zero elements found')

        right_lim += 5
        expectation_data = expectation_data[:right_lim]
        experimental_data = experimental_data[:right_lim]

    if truncate_left:
        for i, ((x1, y1), (x2, y2)) in enumerate(zip(experimental_data, expectation_data)):
            if y1 >= 1 or y2 >= 1:
                left_lim = i
                break
        else:
            raise ValueError('no non-zero elements found')

        left_lim = max(0, left_lim - 5)
        expectation_data = expectation_data[left_lim:]
        experimental_data = experimental_data[left_lim:]

    if use_log:
        expectation_data = [(x, log(y, base=10)) for (x, y) in expectation_data]
        experimental_data = [(x, log(y, base=10)) for (x, y) in experimental_data]

    expected_label = 'expected values' if dist == 'uniform' else 'uniform distribution'
    experimental_label = 'experimental values' if dist == 'uniform' else 'thermodynamic distribution'
    plot_expectation = list_plot(expectation_data, color='green', size=size, legend_label=expected_label)
    plot_experimental = list_plot(experimental_data, color='red', size=size, legend_label=experimental_label)

    max_y = max(itertools.chain(expectation_data, experimental_data), key=lambda xy: xy[1])[1]
    min_y = 0

    my_plot = plot_expectation + plot_experimental
    my_plot.set_axes_range(ymin=min_y, ymax=max_y)
    my_plot.axes_labels([xlabel, 'log of frequency' if use_log else 'frequency'])

    labels = my_plot._extra_kwds.get('axes_labels')
    displace = -.18
    if labels and labels[0]:
        displace -= .005 * my_plot.fontsize()

    my_plot.set_legend_options(
        bbox_to_anchor=(0., displace, 1., .102),
        loc=3,
        mode='expand',
        borderaxespad=0.,
        markerscale=8.0/size,
    )

    return my_plot


def make_plots(base_name, stat_prefix, args, size, xlabel, dist, truncate_left=False, truncate_right=False):
    source_name = os.path.join('data', 'processed_plot_data', 'plotData_' + stat_prefix + base_name)
    with open(source_name, 'r') as f:
        f.readline() # header
        experimental_data, expectation_data = zip(*(map(float, line.split()) for line in f))
    print('experimental data:\n{}\nexpectation:\n{}\n'.format(experimental_data[:10], expectation_data[:10]))
    
    experimental_data = list(experimental_data)
    expectation_data = list(expectation_data)

    reg_plot = construct_plot(expectation_data, experimental_data, size, xlabel, use_log=False, dist=dist, truncate_left=truncate_left, truncate_right=truncate_right)
    log_plot = construct_plot(expectation_data, experimental_data, size, xlabel, use_log=True, dist=dist, truncate_left=truncate_left, truncate_right=truncate_right)

    char_display_name = {
        'cd_sums_' : 'contact distance sums',
        'num_leaves_' : 'number of leaves',
        'root_degree_' : 'root degree',
        'height_' : 'height',
        'ladder_distance_' : 'ladder distance',
        'cd_averages_' : 'contact distance averages',
    }

    no_expected_data = (stat_prefix == 'ladder_distance_') # don't have expectations for ladder distance so this plot is of two sets of experimental data

    for (my_plot, log_prefix) in ((reg_plot, ''), (log_plot, 'log_')):
        plot_name = stat_prefix + log_prefix + base_name[:-4] + '.png'
        plot_path = os.path.join('plots', plot_name)
        print('saving {} to: {}'.format(stat_prefix[:-1], plot_path))
        dist_to_use = 'uniform' if dist == 'uniform' else 'thermodynamic'
        other_dist = 'thermodynamic' if dist_to_use == 'uniform' else 'uniform'
        if no_expected_data:
            plot_title = 'frequencies of {} under \n{} distribution vs. {} distribution'.format(
                char_display_name[stat_prefix], dist_to_use, other_dist)
        else:
            plot_title = 'frequencies of {} under {} distribution\nvs. expected frequencies under uniform distribution'.format(
                char_display_name[stat_prefix], dist_to_use)
        params = 'with n={:,}, initial mixing time of {:,}, sample\ninterval of {:,}, and {:,} samples'.format(
            args.n, args.mixing_time, args.sample_interval, args.num_samples)
        my_plot.save(plot_path, title='\n '.join([plot_title, params]) + '\n\n')


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

    base_name = 'n={}_dist={}_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(args.n, distribution, args.mixing_time, args.sample_interval, args.num_samples)

    if distribution == 'nntm':
        make_plots(base_name, 'ladder_distance_', args, 5, 'ladder distance', distribution, truncate_right=True)

    make_plots(base_name, 'height_', args, 5, 'height', distribution, truncate_right=True)
    make_plots(base_name, 'cd_sums_', args, 5, 'contact distance', distribution)
    make_plots(base_name, 'num_leaves_', args, 5, 'number of leaves', distribution, truncate_right=True, truncate_left=True)
    make_plots(base_name, 'root_degree_', args, 20, 'root degree', distribution, truncate_right=True)
