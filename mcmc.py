# The program will know the number of edges n in plane trees that we are runnining
# Will take in an initial dyckWord
# An energy function E(T) is set
# The markov chain will run t-steps
# While the number of steps is still less than t:
    # From the current plane tree, pick a candidate new dyckword to move to
    # Compute chanceOfMove=min (1, exp(-new_word)/exp(-curr_word)) <---- this is MCMC
    # with prob chanceOfMove, set newDyckWord to be our current dyckWord, otherwise keep current dyckWord as my state

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
        probability = min(1, np.exp(-new_energy)/np.exp(-curr_energy)) # MCMC part
        ran = random.random()
        if ran <= probability: # leave new dyck_word as current state
            return
        else: # revert to old dyck_word
            curr_word[pos1], curr_word[pos0] = curr_word[pos0], curr_word[pos1] # swap back


# energy function rewritten by Anna to be faster and avoid recursion depth issues
def calculate_energy(word):
    root_deg = 0
    num_leaves = 0
    int_nodes = 0
    num_edges = len(word) / 2
    curr_depth = 0
    candidate_int_node_depths = []
    prev_letter = 0 # careful here if you add more details to energy function
    
    for letter in word:
        if candidate_int_node_depths != [] and candidate_int_node_depths[len(candidate_int_node_depths) - 1] == curr_depth:
                    candidate = candidate_int_node_depths.pop()
        if letter == 1: # open of arc
            if curr_depth==0:
                root_deg+=1
            if prev_letter == 1: # two opens in a row
                candidate_int_node_depths.append(curr_depth)
            curr_depth+=1
        else: # close
            if prev_letter == 1: # leaf
                num_leaves+=1
            else: # two closes in a row
                if candidate == curr_depth: # internal node
                    int_nodes+=1
            curr_depth-=1
        prev_letter = letter
        candidate=-1

    return (-.4*root_deg + 2.3*num_leaves +1.3*int_nodes) - .1*num_edges


def contact_distances(word): # For any dyck word, returns the contact distances of the corresponding matching
    # formatted as a list l of length 2n, where l[i] is the number of pairs of cd i
    l = len(word)
    cds = [0 for i in range(l)]
    unpaired = [] # stack of positions of unpaired ('s
    for pos, letter in enumerate(word):
        if letter == 1: # start of an arc
            unpaired.append(pos)
        else: # end of arc
            begin = unpaired.pop()
            dist = pos-begin-1
            cds[dist] += 1
    return cds


def num_leaves(word):
    leaves = 0
    for char1, char2 in zip(word[:-1], word[1:]):
        if (char1 == 1) and (char2 == 0):
            leaves += 1
    return leaves


def root_degree(word):
    degree = 0
    curr_depth = 0
    for char in word:
        if char == 1:
            if curr_depth == 0:
                degree += 1
            curr_depth += 1
        else:
            curr_depth -= 1
    return degree


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


def my_project(start_word, mixing_time, sample_interval, num_samples, distribution):
    """
    start_word: word to start with
    mixing_time: t
    sample_interval: collect every x-amount of steps
    num_samples: number of samples that I want
    distribution: one of 'uniform' or 'nntm', passed to combined_move
    """
    assert distribution in ('uniform', 'nntm')

    print('running my_project with mixing_time={}, sample_interval={}, num_samples={}'.format(
        mixing_time,
        sample_interval,
        num_samples))

    n = len(start_word) / 2
    # samples = [] # an empty list that will append the samples
    curr_word = start_word
    # curr_energy = get_arcs_init_energy(len(start_word)/2)
    for i in range(mixing_time):  # I need my movingWithProb to run mixing_time amount of times (while loop)
        combined_move(curr_word, distribution)
    samp_count = 0
    step_count = 0
    checkpoint = int(num_samples / 10)

    cd_sums = [0] * (2 * n)
    num_leaves_frequency = [0] * n
    root_degree_frequency = [0] * n
    height_frequency = [0] * n

    ### one entry is filled in each time a sample is collected
    num_leaves_values = np.zeros(num_samples, dtype='uint8')
    root_degree_values = np.zeros(num_samples, dtype='uint8')
    height_values = np.zeros(num_samples, dtype='uint8')

    while samp_count < num_samples:  # I need my program to stop after I have collected num_samples amount of samples
        if step_count == sample_interval: # after sample_interval amount of steps, append the curr_word to my list 'samples' (for loop)
            # samples.append(list(curr_word))
            update_frequencies(curr_word, num_leaves_frequency, root_degree_frequency, height_frequency) ### calculate updated statistics for curr_word
            update_samples(curr_word, samp_count, num_leaves_values, root_degree_values, height_values)
            cds = contact_distances(curr_word)
            cd_sums = map(add, cd_sums, cds)

            step_count = 0
            samp_count += 1
            if samp_count % checkpoint == 0:
                print('collected {} of {} samples'.format(samp_count, num_samples))
        else:
            combined_move(curr_word, distribution)
            step_count += 1
    return {
        'frequencies' : {
            'cd_sums' : cd_sums,
            'num_leaves' : num_leaves_frequency,
            'root_degree' : root_degree_frequency,
            'height' : height_frequency
        },
        'samples' : {
            'num_leaves' : num_leaves_values,
            'root_degree' : root_degree_values,
            'height' : height_values
        }}


def update_frequencies(word, num_leaves_frequency, root_degree_frequency, height_frequency):
    leaves = num_leaves(word)
    num_leaves_frequency[leaves] += 1

    degree = root_degree(word)
    root_degree_frequency[degree] += 1

    tree_height = height(word)
    height_frequency[tree_height] += 1


def update_samples(word, i, num_leaves_values, root_degree_values, height_values):
    leaves = num_leaves(word)
    num_leaves_values[i] = leaves

    degree = root_degree(word)
    root_degree_values[i] = degree

    tree_height = height(word)
    height_values[i] = tree_height


def write_to_file(data, base_name, prefix, *dirs):
        filename = os.path.join(os.path.join(*dirs), prefix + base_name) ### don't judge me
        print('saving {} to: {}'.format(prefix[:-1], filename))
        with open(filename, 'w') as f:
            f.writelines('{}\n'.format(i) for i in data)


def write_frequencies_to_file(data, base_name, prefix):
    return write_to_file(data, base_name, prefix, 'data', 'by_frequency')


def write_samples_to_file(data, base_name, prefix):
    return write_to_file(data, base_name, prefix, 'data', 'by_sample')


if __name__ == '__main__':
    np.seterr('raise')
    random.seed(1)

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

    start_time = time.time()

    start_word = [1]*args.n + [0]*args.n

    results = my_project(start_word, args.mixing_time, args.sample_interval, args.num_samples, distribution)
    cd_sums = results['frequencies']['cd_sums']
    num_leaves_frequency = results['frequencies']['num_leaves']
    root_degree_frequency = results['frequencies']['root_degree']
    height_frequency = results['frequencies']['height']

    ### np ndarrays
    num_leaves_values = results['samples']['num_leaves']
    root_degree_values = results['samples']['root_degree']
    height_values = results['samples']['height']


    end_time = time.time()
    print('Elapsed time was {:.0f} seconds.'.format(end_time - start_time))

    cd_sums = list(cd_sums)

    base_name = 'n={}_dist={}_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(args.n, distribution, args.mixing_time, args.sample_interval, args.num_samples)
    write_frequencies_to_file(cd_sums, base_name, 'cd_sums_')
    write_frequencies_to_file(num_leaves_frequency, base_name, 'num_leaves_')
    write_frequencies_to_file(root_degree_frequency, base_name, 'root_degree_')
    write_frequencies_to_file(height_frequency, base_name, 'height_')

    for sample_values_array, prefix in (
        (num_leaves_values, 'num_leaves_'),
        (root_degree_values, 'root_degree_'),
        (height_values, 'height_')):

        filename = os.path.join('data', 'by_sample', prefix + base_name)
        np.savetxt(filename, sample_values_array, fmt='%d')
