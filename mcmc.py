# The program will know the number of edges n in plane trees that we are runnining
# Will take in an initial dyckWord
# An energy function E(T) is set
# The markov chain will run t-steps
# While the number of steps is still less than t:
    # From the current plane tree, pick a candidate new dyckword to move to
    # Compute chanceOfMove=min (1, exp(-newWord)/exp(-currWord)) <---- this is MCMC
    # with prob chanceOfMove, set newDyckWord to be our current dyckWord, otherwise keep current dyckWord as my state

import matplotlib.pyplot as plt
import numpy as np

from functools import lru_cache
import math
import pdb
import random
import sys
import time
import timeit
from operator import add

opt_num = 1

def isValid(dyckWord): #this function checks for valid dyckWords/ballot sequences
    counter1=0
    counter0=0
    for letter in dyckWord:
        if letter ==1:
            counter1+=1
        else:
            counter0+=1
        if counter0> counter1:
            return False 
    return True

def findPos(i, j, dyckWord):
    """
    finds position of ith 1 and jth 0
    """
    counter1 = 0 # keep track of how many 1's have been found
    counter0 = 0 # keep track of how many 0's have been found
    pos1 = -1
    pos0 = -1
    for index, letter in enumerate(dyckWord):
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

def combinedMove(currWord):
    currWordEnergy = fasterEnergyFunction(currWord) # calculates energy before modifying currWord

    length=len(currWord) #pick random i,j between 1 and length
    i=random.randrange(1,length/2)
    j=random.randrange(1,length/2) # TODO: why is this only [1, length/2) ?
    pos1, pos0 = findPos(i, j, currWord)
    currWord[pos1], currWord[pos0] = currWord[pos0], currWord[pos1] # swap values; currWord is now newWord

    if not isValid(currWord):
        currWord[pos1], currWord[pos0] = currWord[pos0], currWord[pos1] # no move to make; swap back
        return

    newWordEnergy = fasterEnergyFunction(currWord) # calculates energy
    probability= min(1, np.exp(-newWordEnergy)/np.exp(-currWordEnergy)) #MCMC Part!!!!
    ran=random.random()
    if ran <= probability:  # leave new dyckword as current state
        return
    else: # revert to old dyckWord
        currWord[pos1], currWord[pos0] = currWord[pos0], currWord[pos1] # swap back

def myProject(startWord, mixingTimeT, sampleInterval, numOfSamples):
    """
    startWord: word to start with
    mixingTimeT: t
    sampleInterval: collect every x-amount of steps
    numOfSamples: number of samples that I want
    """
    samples=[]#an empty list that will append the samples
    currWord=startWord
    for i in range(mixingTimeT):  #I need my movingWithProb to run mixingTimeT amount of times (while loop)
        # newWord=markovMove(currWord)
        # currWord=movingWithProb(currWord,newWord)
        combinedMove(currWord)
    sampCount=0
    stepCount=0
    while sampCount < numOfSamples:  #I need my program to stop after I have collected numOfSamples amount of samples
        if stepCount==sampleInterval: #after sampleInterval amount of steps, append the currWord to my list 'samples' (for loop)
            samples.append(list(currWord))
            stepCount=0
            sampCount+=1
        else:
            # newWord=markovMove(currWord)
            # currWord=movingWithProb(currWord,newWord)
            combinedMove(currWord)
            stepCount+=1
    return samples #return samples

# energy function rewritten by Anna to be faster and avoid recursion depth issues
def fasterEnergyFunction(word, cache={}):
    tup = tuple(word)
    if tup in cache:
        return cache[tup]

    root_deg = 0
    num_leaves = 0
    int_nodes = 0
    num_edges = len(word) / 2
    curr_depth = 0
    candidate_int_node_depths = []
    prev_letter = 0 #caereful here if you add more details to energy function
    
    for letter in word:
        if candidate_int_node_depths != [] and candidate_int_node_depths[len(candidate_int_node_depths) - 1] == curr_depth:
                    candidate = candidate_int_node_depths.pop()
        if letter == 1: #open of arc
            if curr_depth==0:
                root_deg+=1
            if prev_letter == 1: #two opens in a row
                candidate_int_node_depths.append(curr_depth)
            curr_depth+=1
        else: #close
            if prev_letter == 1: #leaf
                num_leaves+=1
            else: # two closes in a row
                if candidate == curr_depth: #internal node
                    int_nodes+=1
            curr_depth-=1
        prev_letter = letter
        candidate=-1

    answer = (-.4*root_deg +(2.3*num_leaves) +1.3*int_nodes) - .1*num_edges
    cache[tup] = answer
    return answer


def contactDistances(word): #For any dyck word, returns the contact distances of the corresponding matching
    # formatted as a list l of length 2n, where l[i] is the number of pairs of cd i
    l = len(word)
    cds = [0 for i in range(l)]
    unpaired = [] #stack of positions of unpaired ('s
    for pos, letter in enumerate(word):
        if letter == 1: #start of an arc
            unpaired.append(pos)
        else: #end of arc
            begin = unpaired.pop()
            dist = pos-begin-1
            cds[dist]+=1
    return cds

start_time = time.time()

n=100
mixingTimeT, sampleInterval, numOfSamples = 50000, 2000, 100
# outPutSamples= myProject(startWord, 1000, 1000, 6)
startWord = [1]*n + [0]*n
outPutSamples= myProject(startWord, mixingTimeT, sampleInterval, numOfSamples)

end_time = time.time()
print("Elapsed time was {:.0f} seconds.".format(end_time - start_time))

#compute contact distances
cd_sums = [0 for i in range(2*n)]
for sample in outPutSamples:
    cds = contactDistances(sample)
    cd_sums = map(add, cd_sums, cds)

cd_sums = list(cd_sums)

# p1 = list_plot(cd_sums[4:], color='red', size=5)
plt.scatter(range(1, len(cd_sums)-3), cd_sums[4:])
plt.savefig("cds_MCMC_thermo_{}n_{}ini_{}int_afterOpt{}.png".format(n, mixingTimeT, sampleInterval, opt_num))

# out1 = open('cds_ MCMC_thermo_1000n_50000ini_2000int.txt', 'w')
# for d, s in enumerate(cd_sums):
#     out1.write(str(d) + "\t" + str(s) + "\n")
# out1.close()


# out_file = open('MCMC_thermo_1000n_50000ini_2000int.txt', 'w')
# for s in outPutSamples:
#     out_file.write(str(s)+"\n")
    
# out_file.close()

# def S(n,d): # cd for uniform distribution
#     if d%2==1:
#         return 0
#     return 1 / (d/2 +1) * binomial(d, d/2) * binomial(2*n - d -1, n - d/2 -1)

# r = 100 / catalan_number(1000)
# p2 = list_plot([r*S(1000, d) for d in range(4, 2000, 1)], color='green', size=5)
# show(p1+p2)
# print [r*S(1000, d).n() for d in range(0, 20, 1)]

### for testing contactDistances
# test_word = DyckWord([1,1,1,0,1,0,0,1,0,0])
# contactDistances(test_word)
