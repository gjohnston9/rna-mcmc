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


def find_pos(i, j, dyck_word):
    """
    finds position of ith 1 and jth 0
    """
    counter1 = 0 # keep track of how many 1's have been found
    counter0 = 0 # keep track of how many 0's have been found
    pos1 = -1
    pos0 = -1
    for index, letter in enumerate(dyck_word):
        if letter == 1:
            counter1 += 1
            if counter1 == i:
                pos1 = index # we found the ith 1
        else:
            counter0 += 1
            if counter0 == j:
                pos0 = index # we found the jth 0
        if -1 not in [pos1, pos0]:
            break

    assert -1 not in [pos1, pos0]
    return pos1, pos0


def combined_move(curr_word, move_type, distribution):
    length = len(curr_word) # pick random i,j between 1 and length

    if move_type == 'by_count':
        i = random.randrange(1, length/2 + 1)
        j = random.randrange(1, length/2 + 1)
        pos1, pos0 = findPos(i, j, curr_word)
    elif move_type == 'by_position':
        pos1 = random.randrange(0, length)
        pos0 = random.randrange(0, length)
        if curr_word[pos1] == curr_word[pos0]: # self-loop
            return

    if distribution == 'nntm':
        curr_energy = calculate_energy(curr_word) # calculates energy before modifying curr_word

    curr_word[pos1], curr_word[pos0] = curr_word[pos0], curr_word[pos1] # swap values; curr_word is now new_word

    if not is_valid(curr_word, end=max(pos1, pos0)):
        curr_word[pos1], curr_word[pos0] = curr_word[pos0], curr_word[pos1] # no move to make; swap back
        return

    if distribution == 'nntm': # if distribution is uniform, we make the above move 100% of the time (never swap back)
        new_energy = calculate_energy(curr_word) # calculates energy
        probability = min(1, np.exp(-new_energy)/np.exp(-curr_energy)) # MCMC part
        ran = random.random()
        if ran <= probability:  # leave new dyck_word as current state
            return
        else: # revert to old dyck_word
            curr_word[pos1], curr_word[pos0] = curr_word[pos0], curr_word[pos1] # swap back


def my_project(start_word, mixing_time, sample_interval, num_samples, move_type, distribution):
    """
    start_word: word to start with
    mixing_time: t
    sample_interval: collect every x-amount of steps
    num_samples: number of samples that I want
    move_type: one of 'by_count' or 'by_position', passed to combined_move
    distribution: one of 'uniform' or 'nntm', passed to combined_move
    """
    assert move_type in ('by_count', 'by_position')
    assert distribution in ('uniform', 'nntm')

    samples = [] # an empty list that will append the samples
    curr_word = start_word
    for i in range(mixing_time):  # I need my movingWithProb to run mixing_time amount of times (while loop)
        combined_move(curr_word, move_type, distribution)
    samp_count = 0
    step_count = 0
    while samp_count < num_samples:  # I need my program to stop after I have collected num_samples amount of samples
        if step_count == sample_interval: # after sample_interval amount of steps, append the curr_word to my list 'samples' (for loop)
            samples.append(list(curr_word))
            step_count = 0
            samp_count += 1
        else:
            combined_move(curr_word, move_type, distribution)
            step_count += 1
    return samples # return samples


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

    return (-.4*root_deg +(2.3*num_leaves) +1.3*int_nodes) - .1*num_edges


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


if __name__ == '__main__':
    np.seterr('raise')

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

    start_time = time.time()

    start_word = [1] * args.n + [0] * args.n
    output_samples = my_project(start_word, args.mixing_time, args.sample_interval, args.num_samples, move_type, distribution)

    end_time = time.time()
    print('Elapsed time was {:.0f} seconds.'.format(end_time - start_time))

    # compute contact distances
    cd_sums = [0 for i in range(2 * args.n)]
    for sample in output_samples:
        cds = contact_distances(sample)
        cd_sums = map(add, cd_sums, cds)

    cd_sums = list(cd_sums)

    cd_sums_name = 'sums_n={}_moveType={}_dist={}_mixingTime={}_sample_interval={}_num_samples={}.txt'.format(args.n, move_type, distribution, args.mixing_time, args.sample_interval, args.num_samples)
    cd_sums_name = os.path.join('data', cd_sums_name)
    print('saving cd_sums to: {}'.format(cd_sums_name))

    with open(cd_sums_name, 'w') as f:
        f.writelines('{}\n'.format(i) for i in cd_sums)

    # ### for testing contact_distances
    # # test_word = DyckWord([1,1,1,0,1,0,0,1,0,0])
    # # contact_distances(test_word)
