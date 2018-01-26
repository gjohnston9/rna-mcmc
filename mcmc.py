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
from collections import namedtuple
from operator import add

arc_tuple = namedtuple('Arc', ['partner', 'first_child', 'next_sibling', 'parent'])

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

    return

    # ### posA and posB hold the INDICES IN CURR_WORD (0-based) of swapped chars
    # endpointA = posA + 1
    # endpointB = posB + 1
    # ### now they refer to endpoints

    # if (endpointA + endpointB) % 2 == 1: ### one is odd, one is even; they may refer to endpoints of the same arc (in orig curr_word)
    #     if endpointB % 2 == 1:
    #         endpointA, endpointB = endpointB, endpointA
    #     ### endpointA is odd, endpointB is even
    #     indexA = arc_to_index(endpointA)
    #     arcA = curr_arcs[indexA]
    #     if arcA.partner == endpointB:
    #         assert arcA.parent != None
    #         endpointB = arcA.parent
    #         indexB = arc_to_index(endpointB)
    #         ### now endpointA and endpointB refer to odd endpoints of two unobstructed arcs,
    #         ### and indexA and indexB refer to the indices of these arcs in curr_arcs
    #     else:
    #         ### endpointB refers to the even endpoint of some other arc. Time to find it
    #         for i, arc in enumerate(curr_arcs):
    #             if arc.partner == endpointB:
    #                 indexB = i
    #                 endpointB = index_to_arc(indexB)
    #                 break
    #         else:
    #             raise ValueError('could not find arc with even endpoint {}'.format(posB))
    # elif endpointA % 2 == 0: ### they are both even
    #     foundA, foundB = False, False
    #     for i, arc in enumerate(curr_arcs):
    #         if arc.partner == endpointA:
    #             indexA = i
    #             foundA = True
    #         elif arc.partner == endpointB:
    #             indexB = i
    #             foundB = True
    #         if foundA and foundB:
    #             break
    #     else:
    #         err = ''
    #         if not foundA: err += 'could not find arc with endpointA ({})\n'.format(endpointA)
    #         if not foundB: err += 'could not find arc with endpointB ({})\n'.format(endpointB)
    #         raise ValueError(err)
    # else: ### they are both odd
    #     indexA = arc_to_index(endpointA)
    #     indexB = arc_to_index(endpointB)

    curr_arcs = curr_state
    ### pick the odd endpoints of two arcs to potentially interchange
    length = len(curr_arcs)
    index0 = random.randrange(0, length)
    index1 = random.randrange(0, length)
    # index0 = indexA
    # index1 = indexB
    index0, index1 = min(index0, index1), max(index0, index1)
    arc0_endpoint = index_to_arc(index0)
    arc1_endpoint = index_to_arc(index1)
    arc0 = curr_arcs[index0]
    arc1 = curr_arcs[index1]

    siblings = False
    arc0_child_of_arc1 = False
    arc1_child_of_arc0 = False

    if index0 == index1:
        return None

    # print('index0: {}\nindex1: {}\narc0_endpoint: {}\narc1_endpoint: {}\n\n'.format(index0, index1, arc0_endpoint, arc1_endpoint))

    if arc0.parent == arc1.parent:
        siblings = True
    if arc0.parent == arc1_endpoint:
        arc0_child_of_arc1 = True
    if arc1.parent == arc0_endpoint:
        arc1_child_of_arc0 = True

    s = sum([siblings, arc0_child_of_arc1, arc1_child_of_arc0])
    if s == 0:
        return None
    elif s != 1:
        pdb.set_trace()

    try:
        if siblings:
            if arc0_endpoint < arc0.partner:
                assert arc1.partner > arc1_endpoint
                arc0_new_first_child = arc0.first_child or arc1_endpoint
                curr_arcs[index0] = arc_tuple(
                    partner = arc1.partner,
                    first_child = arc0_new_first_child,
                    next_sibling = arc1.next_sibling,
                    parent = arc0.parent)
                curr_arcs[index1] = arc_tuple(
                    partner = arc0.partner,
                    first_child = arc0.next_sibling if arc0.next_sibling < arc1_endpoint else None,
                    next_sibling = arc1.first_child,
                    parent = arc0_endpoint)

                start_points_and_end_points_and_new_parents = (
                    (arc0.next_sibling, arc1_endpoint, arc1_endpoint),
                    (arc1.first_child, None, arc0_endpoint))

                gains_sibling_start_endpoint = arc0.first_child ### start from here and follow sibling path
                gains_sibling_new_sibling_endpoint = arc1_endpoint
                loses_sibling_start_endpoint = arc0.next_sibling
                loses_sibling_end_endpoint = arc1_endpoint

            else: ### arc0_endpoint > arc0.partner
                assert arc1.partner < arc1_endpoint
                arc1_new_first_child = arc0.first_child or arc0_endpoint
                curr_arcs[index0] = arc_tuple(
                    partner = arc1.partner,
                    first_child = arc0.next_sibling if arc0.next_sibling < arc1_endpoint else None,
                    next_sibling = arc1.first_child,
                    parent = arc1_endpoint)
                curr_arcs[index1] = arc_tuple(
                    partner = arc0.partner,
                    first_child = arc1_new_first_child,
                    next_sibling = arc1.next_sibling,
                    parent = arc1.parent)

                ### have to update parent arc's first child if it was originally arc0
                parent_arc_endpoint = arc0.parent
                parent_index = arc_to_index(parent_arc_endpoint)
                parent_arc = curr_arcs[parent_index]
                if parent_arc.first_child == arc0_endpoint:
                    curr_arcs[parent_index] = arc_tuple(
                        partner = parent_arc.partner,
                        first_child = arc1_endpoint,
                        next_sibling = parent_arc.next_sibling,
                        parent = parent_arc.parent)
                else:
                    ### otherwise, original arc0's previous sibling's next_sibling has changed
                    curr_arc_endpoint = parent_arc.first_child
                    while curr_arc_endpoint != arc0_endpoint:
                        curr_index = arc_to_index(curr_arc_endpoint)
                        curr_arc = curr_arcs[curr_index]
                        curr_arc_endpoint = curr_arc.next_sibling
                    ### curr_index and curr_arc now refer to orig arc0's prev sibling
                    assert curr_arc.next_sibling == arc0_endpoint
                    curr_arcs[curr_index] = arc_tuple(
                        partner = curr_arc.partner,
                        first_child = curr_arc.first_child,
                        next_sibling = arc1_endpoint, ### update next_sibling
                        parent = curr_arc.parent)

                start_points_and_end_points_and_new_parents = (
                    (arc0.first_child, None, arc1_endpoint),
                    (arc0.next_sibling, arc1_endpoint, arc0_endpoint))

                gains_sibling_start_endpoint = arc0.first_child ### start from here and follow sibling path
                gains_sibling_new_sibling_endpoint = arc0_endpoint
                loses_sibling_start_endpoint = arc0.next_sibling
                loses_sibling_end_endpoint = arc1_endpoint

        elif arc0_child_of_arc1:
            assert arc0.partner > arc0_endpoint
            assert arc1.partner < arc1_endpoint
            arc0_new_first_child = arc1.first_child if (arc1.first_child < arc0_endpoint) else None
            arc0_new_sibling = arc0.first_child or arc1_endpoint
            curr_arcs[index0] = arc_tuple(
                partner = arc1.partner,
                first_child = arc0_new_first_child,
                next_sibling = arc0_new_sibling,
                parent = arc1.parent)
            curr_arcs[index1] = arc_tuple(
                partner = arc0.partner,
                first_child = arc0.next_sibling,
                next_sibling = arc1.next_sibling,
                parent = arc1.parent)

            ### have to update arc1's parent arc's first child if it was originally arc1
            parent_arc_endpoint = arc1.parent
            parent_index = arc_to_index(parent_arc_endpoint)
            parent_arc = curr_arcs[parent_index]
            if parent_arc.first_child == arc1_endpoint:
                curr_arcs[parent_index] = arc_tuple(
                    partner = parent_arc.partner,
                    first_child = arc0_endpoint,
                    next_sibling = parent_arc.next_sibling,
                    parent = parent_arc.parent)
            else:
                ### otherwise, original arc1's previous sibling's next_sibling has changed
                curr_arc_endpoint = parent_arc.first_child
                while curr_arc_endpoint != arc1_endpoint:
                    curr_index = arc_to_index(curr_arc_endpoint)
                    curr_arc = curr_arcs[curr_index]
                    curr_arc_endpoint = curr_arc.next_sibling
                ### curr_index and curr_arc now refer to orig arc1's prev sibling
                assert curr_arc.next_sibling == arc1_endpoint
                curr_arcs[curr_index] = arc_tuple(
                    partner = curr_arc.partner,
                    first_child = curr_arc.first_child,
                    next_sibling = arc0_endpoint, ### update next_sibling
                    parent = curr_arc.parent)

            start_points_and_end_points_and_new_parents = (
                (arc1.first_child, arc0_endpoint, arc0_endpoint),
                (arc0.first_child, None, arc1.parent))

            gains_sibling_start_endpoint = arc0.first_child
            gains_sibling_new_sibling_endpoint = arc1_endpoint
            loses_sibling_start_endpoint = arc1.first_child
            loses_sibling_end_endpoint = arc0_endpoint

        else: ### arc1_child_of_arc0
            assert arc0.partner > arc0_endpoint
            assert arc1.partner < arc1_endpoint
            arc0_new_first_child = arc0.first_child if (arc0.first_child < arc1_endpoint) else None
            arc0_new_sibling = arc1.first_child or arc1_endpoint
            curr_arcs[index0] = arc_tuple(
                partner = arc1.partner,
                first_child = arc0_new_first_child,
                next_sibling = arc0_new_sibling,
                parent = arc0.parent)
            curr_arcs[index1] = arc_tuple(
                partner = arc0.partner,
                first_child = arc1.next_sibling,
                next_sibling = arc0.next_sibling,
                parent = arc0.parent)

            start_points_and_end_points_and_new_parents = (
                (arc1.first_child, None, arc0.parent),
                (arc1.next_sibling, None, arc1_endpoint))

            gains_sibling_start_endpoint = arc1.first_child
            gains_sibling_new_sibling_endpoint = arc1_endpoint
            loses_sibling_start_endpoint = arc0.first_child
            loses_sibling_end_endpoint = arc1_endpoint

        for start_point, end_point, new_parent in start_points_and_end_points_and_new_parents:
            curr_arc_endpoint = start_point
            while (curr_arc_endpoint is not None) and (curr_arc_endpoint != end_point):
                index = arc_to_index(curr_arc_endpoint)
                curr_arc = curr_arcs[index]
                curr_arcs[index] = arc_tuple(
                    partner = curr_arc.partner,
                    first_child = curr_arc.first_child,
                    next_sibling = curr_arc.next_sibling,
                    parent = new_parent) ### update parent
                curr_arc_endpoint = curr_arc.next_sibling

        curr_arc_endpoint = gains_sibling_start_endpoint
        curr_arc = None
        index = None
        while curr_arc_endpoint is not None:
            index = arc_to_index(curr_arc_endpoint)
            curr_arc = curr_arcs[index]
            curr_arc_endpoint = curr_arc.next_sibling
        if curr_arc is not None:
            curr_arcs[index] = arc_tuple(
                partner = curr_arc.partner,
                first_child = curr_arc.first_child,
                next_sibling = gains_sibling_new_sibling_endpoint, ### update next_sibling
                parent = curr_arc.parent)

        curr_arc_endpoint = loses_sibling_start_endpoint
        curr_arc = None
        index = None
        while (curr_arc_endpoint is not None) and (curr_arc_endpoint != loses_sibling_end_endpoint):
            index = arc_to_index(curr_arc_endpoint)
            curr_arc = curr_arcs[index]
            curr_arc_endpoint = curr_arc.next_sibling
        if curr_arc is not None:
            curr_arcs[index] = arc_tuple(
                partner = curr_arc.partner,
                first_child = curr_arc.first_child,
                next_sibling = None,
                parent = curr_arc.parent)            

    except AssertionError:
        pdb.set_trace()

    return 1 ### success


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

    samples = [] # an empty list that will append the samples
    curr_word = start_word
    # curr_energy = get_arcs_init_energy(len(start_word)/2)
    for i in range(mixing_time):  # I need my movingWithProb to run mixing_time amount of times (while loop)
        combined_move(curr_word, distribution)
    samp_count = 0
    step_count = 0
    while samp_count < num_samples:  # I need my program to stop after I have collected num_samples amount of samples
        if step_count == sample_interval: # after sample_interval amount of steps, append the curr_word to my list 'samples' (for loop)
            samples.append(list(curr_word))
            step_count = 0
            samp_count += 1
        else:
            combined_move(curr_word, distribution)
            step_count += 1
    return samples


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


def num_leaves(word):
    leaves = 0
    for i, char in enumerate(word[:-1]):
        if (char == 1) and word[i+1] == 0:
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


# def index_to_arc(num):
#     """
#     takes index (0-based) of the arc in our representation with odd endpoint
#     n (starting at 1) and returns n
#     """
#     return 2*num + 1


# def arc_to_index(num):
#     """
#     takes odd arc endpoint and returns index of that arc in our representation
#     """
#     assert num % 2 == 1
#     return (num - 1)/2


# def get_arcs_flat_representation(n):
#     """
#     Returns an alternate representation of [1,0]*n. This representation is a list
#     of n elements, each of which is a list of four elements that describes an arc.
    
#     The list at index i in the returned list corresponds to the arc with
#     odd endpoint 2i + 1, and contains these four elements:
#     [
#         the even endpoint of this arc,
#         the odd endpoint of the first child of this arc, or None if this arc has no children,
#         the odd endpoint of the next sibling of this arc, or None if this arc has no next sibling,
#         the odd endpoint of the parent of this arc, or None if this arc has no parent
#     ]
#     """
#     return [arc_tuple(2*i, None, 2*i + 1, None) for i in range(1, n)] + [arc_tuple(2*n, None, None, None)]


# def get_arcs_stacked_representation(n):
#     """
#     Returns arc representation of [1]*n + [0]*n
#     """
#     arcs = [0]*n ### to be replaced by arc_tuples

#     mid_index = int(n/2)
#     mid_odd_endpoint = index_to_arc(mid_index) ### odd endpoint of middle arc
#     partners = range(2*n, 0, -2)
#     children = range(2*n-1, n, -2) + [None] + range(n-1 if n%2==0 else n, 1, -2)
#     parents = [None] + range(2*n-1, mid_odd_endpoint, -2) + range(mid_odd_endpoint-2, 0, -2)

#     for i, (partner, child, parent) in enumerate(zip(partners, children, parents)):
#         arcs[i] = arc_tuple(partner, child, None, parent)
#     return arcs

# def arcs_to_word(arcs):
#     answer = [0]*(2*len(arcs))
#     for index, arc in enumerate(arcs):
#         a = index_to_arc(index) - 1 ### odd endpoint of arc, using 0-indexing
#         b = arc.partner - 1 ### even endpoint of arc, using 0-indexing
#         answer[min(a,b)] = 1
#     return answer


# def get_arcs_init_energy(n):
#     """
#     Returns the initial energy for the structure returned by get_arcs_representation
#     """
#     root_deg = n
#     num_leaves = n
#     int_nodes = 0
#     num_edges = n
#     return -.4*root_deg + 2.3*num_leaves + 1.3*int_nodes - .1*num_edges


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

    output_samples = my_project(start_word, args.mixing_time, args.sample_interval, args.num_samples, distribution)

    end_time = time.time()
    print('Elapsed time was {:.0f} seconds.'.format(end_time - start_time))

    # compute contact distances and num_leaves
    cd_sums = [0] * (2 * args.n)
    num_leaves_frequency = [0] * args.n
    root_degree_frequency = [0] * args.n

    for word in output_samples:
        cds = contact_distances(word)
        cd_sums = map(add, cd_sums, cds)

        leaves = num_leaves(word)
        num_leaves_frequency[leaves] += 1

        degree = root_degree(word)
        root_degree_frequency[degree] += 1

    cd_sums = list(cd_sums)

    cd_sums_name = 'cd_sums_n={}_dist={}_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(args.n, distribution, args.mixing_time, args.sample_interval, args.num_samples)
    cd_sums_name = os.path.join('data', cd_sums_name)
    print('saving cd_sums to: {}'.format(cd_sums_name))

    with open(cd_sums_name, 'w') as f:
        f.writelines('{}\n'.format(i) for i in cd_sums)


    num_leaves_name = 'num_leaves_n={}_dist={}_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(args.n, distribution, args.mixing_time, args.sample_interval, args.num_samples)
    num_leaves_name = os.path.join('data', num_leaves_name)
    print('saving num_leaves to: {}'.format(num_leaves_name))

    with open(num_leaves_name, 'w') as f:
        f.writelines('{}\n'.format(i) for i in num_leaves_frequency)



    root_degree_name = 'root_degree_n={}_dist={}_mixingTime={}_sampleInterval={}_numSamples={}.txt'.format(args.n, distribution, args.mixing_time, args.sample_interval, args.num_samples)
    root_degree_name = os.path.join('data', root_degree_name)
    print('saving root_degree to: {}'.format(root_degree_name))

    with open(root_degree_name, 'w') as f:
        f.writelines('{}\n'.format(i) for i in root_degree_frequency)
