# The program will know the number of edges n in plane trees that we are runnining
# Will take in an initial dyckWord
# An energy function E(T) is set
# The markov chain will run t-steps
# While the number of steps is still less than t:
    # From the current plane tree, pick a candidate new dyckword to move to
    # Compute chanceOfMove=min (1, exp(-new_word)/exp(-curr_word)) <---- this is MCMC
    # with prob chanceOfMove, set newDyckWord to be our current dyckWord, otherwise keep current dyckWord as my state

import numpy as np

import argparse
import math
import os
import pdb
import random
import sys
import time
import timeit
from operator import add


def is_valid(dyck_word, end):
    """
    this function checks for valid dyck_words/ballot sequences, ignoring letters past the `end` index
    """
    counter1 = 0
    counter0 = 0
    for letter in dyck_word[:end+1]:
        if letter == 1:
            counter1 += 1
        else:
            counter0 += 1
        if counter0 > counter1:
            return False 
    return True


def swap(arr, posA, posB):
    arr[posA], arr[posB] = arr[posB], arr[posA]


def combined_move(curr_word, distribution):
    length = len(curr_word)
    
    if distribution == 'nntm':
        old_energy = calculate_energy(curr_word)

    posA = random.randrange(0, length)
    posB = random.randrange(0, length)
    if curr_word[posA] == curr_word[posB]: # self-loop
        return

    swap(curr_word, posA, posB) # curr_arcs is now new_word

    if not is_valid(curr_word, end=max(posA, posB)):
        swap(curr_word, posA, posB) # no move to make; swap back
        return

    if distribution == 'nntm':
        new_energy = calculate_energy(curr_word) # calculates energy
        probability = min(1, np.exp(old_energy-new_energy)) # MCMC part
        assert probability > 0
        ran = random.random()
        if ran <= probability: # leave new dyck_word as current state
            return
        else: # revert to old dyck_word
            swap(curr_word, posA, posB) # curr_arcs is now new_word


def calculate_useful_characteristics(word):
    degree = 0
    leaves = 0
    int_nodes = 0
    num_edges = len(word) / 2
    curr_depth = 0
    candidate_int_node_depths = []
    prev_letter = 0 # careful here if you add more details to energy function
    
    for letter in word:
        if candidate_int_node_depths != [] and candidate_int_node_depths[len(candidate_int_node_depths) - 1] == curr_depth:
                    candidate = candidate_int_node_depths.pop()
        if letter == 1: # open of arc
            if curr_depth == 0:
                degree += 1
            if prev_letter == 1: # two opens in a row
                candidate_int_node_depths.append(curr_depth)
            curr_depth+=1
        else: # close
            if prev_letter == 1: # leaf
                leaves+=1
            else: # two closes in a row
                if candidate == curr_depth: # internal node
                    int_nodes+=1
            curr_depth-=1
        prev_letter = letter
        candidate=-1

    return {
        'degree' : degree,
        'leaves' : leaves,
        'int_nodes' : int_nodes,
        'num_edges' : num_edges,
    }


def calculate_energy(word):
    characteristics = calculate_useful_characteristics(word)
    degree = characteristics['degree']
    leaves = characteristics['leaves']
    int_nodes = characteristics['int_nodes']
    num_edges = characteristics['num_edges']
    return (-.4*degree + 2.3*leaves +1.3*int_nodes) - .1*num_edges


def contact_distances(word): # For any dyck word, returns the contact distances of the corresponding matching
    # formatted as a list l of length 2n, where l[i] is the number of pairs of cd i
    l = len(word)
    cds = [0] * l
    unpaired = [] # stack of positions of unpaired ('s
    for pos, letter in enumerate(word):
        if letter == 1: # start of an arc
            unpaired.append(pos)
        else: # end of arc
            begin = unpaired.pop()
            dist = pos-begin-1
            cds[dist] += 1
    return cds


def contact_distances_average(contact_distances):
    total = 0
    for index, val in enumerate(contact_distances):
        total += index * val
    return total  / (1.0 * len(contact_distances))


def avg_branching(word, root_deg, int_nodes, leaves):
    n = len(word) / 2
    k = 0 if root_deg > 1 else 1
    return (n + 1.0 - int_nodes - k) / (n + 2.0 - int_nodes - leaves - k)


def height(word):
    tree_height = 0
    curr_depth = 0
    for char in word:
        if char == 1:
            curr_depth += 1
            tree_height = max(tree_height, curr_depth)
        else:
            curr_depth -= 1
    return tree_height


def ladder_distance(word):
    curr_depth = 0
    max_depth = 0
    for i, char in enumerate(word):
        if char == 1:
            curr_depth += 1
            if curr_depth > max_depth:
                index_of_deepest_vertex = i
                max_depth = curr_depth
        else:
            curr_depth -= 1

    # now, find the vertex farthest from the deepest vertex
    max_length = 0 # length of path from deepes vertex to the farthest vertex that has been found so far
    for start, direction, close_char in ((index_of_deepest_vertex - 1, -1, 1), (index_of_deepest_vertex + 2, 1, 0)):
        curr_depth = 0
        curr_length = 0
        for char in word[start::direction]:
            if char == close_char:
                if curr_depth == 0:
                    curr_length += 1
                    max_length = max(max_length, curr_length)
                else:
                    curr_depth -= 1
                    curr_length -= 1
            else:
                curr_depth += 1
                curr_length += 1
                max_length = max(max_length, curr_length)

    return max_length + 1


def average_ladder_distance(word):
    """
    The reasoning behind this algorithm (from Anna):
    We want to figure out how much each edge contributes to the average
    ladder distance. That is, we want to compute (the number of paths which
    use the edge)/M. For any edge e, let S be the subtree below e and
    including e itself. Observe that any path joining an edge in S to an
    edge outside of S must pass through edge e. Similarly, any path joining
    two edges in S or two edges in the complement of S will not use edge e,
    unless it is a path joining two edges in S that has e as a terminal edge.
    This case gives us an additional | S | - 1 paths that use e.
    Therefore, the number of paths that use e is | S | * (n minus | S |) + | S | - 1.
    """
    n = len(word) / 2
    total_paths = (n*(n-1)) / 2
    # a 1 at position i will be paired with a 0 at position partners[i]
    partners = [0]*len(word)
    stack = []
    for i, char in enumerate(word):
        if char == 1:
            stack.append(i)
        else:
            partners[stack.pop()] = i

    average_ladder_distance_sum = 0
    for i, char in enumerate(word):
        if char == 0:
            continue
        s = (partners[i] - i + 1)/2 # number of edges in subtree rooted at this edge, including this edge
        average_ladder_distance_sum += ((s)*(n - s) + (s - 1)) / (1.0 * total_paths)

    return average_ladder_distance_sum


def my_project(start_word, mixing_time, sample_interval, num_samples, distribution, base_prefix, base_name):
    """
    start_word: word to start with
    mixing_time: t
    sample_interval: collect every x-amount of steps
    num_samples: number of samples that I want
    distribution: one of 'uniform' or 'nntm', passed to combined_move
    """
    assert distribution in ('uniform', 'nntm')

    print('running my_project with mixing_time={0}, sample_interval={1}, num_samples={2}'.format(
        mixing_time,
        sample_interval,
        num_samples))

    n = len(start_word) / 2
    curr_word = start_word
    # curr_energy = get_arcs_init_energy(len(start_word)/2)
    for i in range(1, mixing_time+1):  # I need my movingWithProb to run mixing_time amount of times (while loop)
        if i > 1 and i % (mixing_time / 10) == 0:
            print('finished {0} of {1} mixing steps'.format(i, mixing_time))
        combined_move(curr_word, distribution)
    samp_count = 0
    step_count = 0
    checkpoint = int(num_samples / 10)

    cd_sums = [0] * (2 * n)
    num_leaves_frequency = np.zeros(n, dtype='uint32')
    root_degree_frequency = np.zeros(n, dtype='uint32')
    height_frequency = np.zeros(n, dtype='uint32')
    ladder_distance_frequency = np.zeros(n, dtype='uint32')

    num_leaves_values = []
    root_degree_values = []
    height_values = []
    ladder_distance_values = []
    branching_values = []
    cd_averages_values = []
    average_ladder_distance_values = []

    batch_size = min(int(1e6), int(num_samples / 10))
    assert num_samples % batch_size == 0 ### don't want to discard any samples

    file_open_mode = 'w' # the first time we open files to write data, we use 'w' to erase whatever was already there
    # it is changed to 'a' after the first time opening files

    while samp_count < num_samples:  # I need my program to stop after I have collected num_samples amount of samples
        if step_count == sample_interval: # after sample_interval amount of steps, append the curr_word to my list 'samples' (for loop)
            ### calculate prelimiary characteristics
            characteristics = calculate_useful_characteristics(curr_word)
            int_nodes = characteristics['int_nodes']
            num_edges = characteristics['num_edges']

            ### calculate characteristics that will be saved
            leaves = characteristics['leaves']
            degree = characteristics['degree']
            tree_height = height(curr_word)
            distance = ladder_distance(curr_word)
            branching = avg_branching(curr_word, degree, int_nodes, leaves)
            cds = contact_distances(curr_word)
            cd_averages = contact_distances_average(cds)
            avg_ladder_distance = average_ladder_distance(curr_word)

            ### update frequencies
            num_leaves_frequency[leaves-1] += 1
            root_degree_frequency[degree-1] += 1
            height_frequency[tree_height-1] += 1
            ladder_distance_frequency[distance-1] += 1
            # no frequency for avg branching, cd_averages or avg_ladder_distance since they are floats
            cd_sums = map(add, cd_sums, cds)

            ### update values
            num_leaves_values.append(leaves)
            root_degree_values.append(degree)
            height_values.append(tree_height)
            ladder_distance_values.append(distance)
            branching_values.append(branching)
            cd_averages_values.append(cd_averages)
            average_ladder_distance_values.append(avg_ladder_distance)

            step_count = 0
            samp_count += 1

            if samp_count % batch_size == 0:
                print('saving batch {0} of {1}'.format(samp_count / batch_size, num_samples / batch_size))
                for array, type_prefix in (
                    (num_leaves_values, 'num_leaves_'),
                    (root_degree_values, 'root_degree_'),
                    (height_values, 'height_'),
                    (ladder_distance_values, 'ladder_distance_'),
                    (branching_values, 'avg_branching_'),
                    (cd_averages_values, 'cd_averages_'),
                    (average_ladder_distance_values, 'avg_ladder_distance_')):

                    filename = os.path.join('data', 'by_sample', base_prefix + type_prefix + base_name)
                    print('saving {0} to {1}'.format(type_prefix[:-1], filename))
                    with open(filename, file_open_mode) as f:
                        for sample in array:
                            f.write('{0}\n'.format(sample))
                num_leaves_values = []
                root_degree_values = []
                height_values = []
                ladder_distance_values = []
                branching_values = []
                cd_averages_values = []
                average_ladder_distance_values = []

                file_open_mode = 'a' # for the rest of the runs, append instead of overwriting files
        else:
            combined_move(curr_word, distribution)
            step_count += 1
    return {
        'cd_sums' : cd_sums,
        'num_leaves' : num_leaves_frequency,
        'root_degree' : root_degree_frequency,
        'height' : height_frequency,
        'ladder_distance' : ladder_distance_frequency,
        }


if __name__ == '__main__':
    np.seterr('raise')

    parser = argparse.ArgumentParser()
    parser.add_argument('n', type=int, help='value of n to use')
    parser.add_argument('mixing_time', type=int,
        help='initial mixing time before recording info about perfect matching characteristics')
    parser.add_argument('sample_interval', type=int,
        help='number of moves to make in between recording info about characteristics')
    parser.add_argument('num_samples', type=int, help='number of samples to collect')

    distribution_selection_group = parser.add_mutually_exclusive_group(required=True)
    distribution_selection_group.add_argument('--uniform', action='store_true',
        help='use a uniform distribution when choosing whether to make a move')
    distribution_selection_group.add_argument('--nntm', action='store_true',
        help='use the ratio of energies as predicted by the nearest-neighbor thermodynamic model ' \
        + 'when choosing whether to make a move')

    parser.add_argument('--prefix', type=str, default='', help='prefix to add to the beginning of ' \
        + 'saved data (useful when repeatedly calling this simulation in a script, so as not to overwrite previous results)')
    parser.add_argument('--start_word_source', type=str, default=None, help='path to file containing ' \
        + 'the starting word to be used (the default start word is [1]*n + [0]*n)')

    args = parser.parse_args()

    if args.uniform:
        distribution = 'uniform'
    else:
        distribution = 'nntm'

    if args.start_word_source is not None:
        with open(args.start_word_source, 'r') as f:
            start_word = list(map(int, f.readline().strip()))
    else:
        start_word = [1] * args.n + [0] * args.n

    base_name = 'n={0}_dist={1}_mixingTime={2}_sampleInterval={3}_numSamples={4}.txt'.format(
        args.n, distribution, args.mixing_time, args.sample_interval, args.num_samples)

    start_time = time.time()
    results = my_project(start_word, args.mixing_time, args.sample_interval, args.num_samples, distribution, args.prefix, base_name)
    end_time = time.time()
    print('Elapsed time was {0:.0f} seconds.'.format(end_time - start_time))
    
    ### np ndarrays
    cd_sums = results['cd_sums']
    cd_sums = np.asarray(list(cd_sums))
    num_leaves_frequency = results['num_leaves']
    root_degree_frequency = results['root_degree']
    height_frequency = results['height']
    ladder_distance_frequency = results['ladder_distance']

    for array, prefix in (
        (cd_sums, 'cd_sums_'),
        (num_leaves_frequency, 'num_leaves_'),
        (root_degree_frequency, 'root_degree_'),
        (height_frequency, 'height_'),
        (ladder_distance_frequency, 'ladder_distance_')):

        filename = args.prefix + prefix + base_name
        filepath = os.path.join('data', 'by_frequency', filename)
        print('saving {0} to data/by_frequency/'.format(filepath))
        np.savetxt(filepath, array, fmt='%d')
