#The program will know the number of edges n in plane trees that we are runnining
#Will take in an initial dyckWord
#An energy function E(T) is set
#The markov chain will run t-steps
#While the number of steps is still less than t:
    #From the current plane tree, pick a candidate new dyckword to move to
    #Compute chanceOfMove=min (1, exp(-newWord)/exp(-currWord)) <---- this is MCMC
    #with prob chanceOfMove, set newDyckWord to be our current dyckWord, otherwise keep current dyckWord as my state

#Complete program!!!!
from sage.combinat.dyck_word import *
import random
# ENERGY FUNCTION
def childrenOfRoot(myTree):
    r=len(myTree)
    return r

def childrenOfVertices(myTree): #this function finds the number of children of each vertex and put in a list
    theList=[] 
    if len(myTree)==0: #base case (also accounts for leaves)
        return [0]
    children=len(myTree) #local calculation
    theList.append(children)
    for i in myTree: 
        subTree=i
        theList.extend(childrenOfVertices(subTree)) #this is where recursion is happening 
    return theList
def returnLeaves(theList):
    d0=0  
    theList=childrenOfVertices(theList)
    theList.pop(0)
    for item in theList:
        if item==0:
            d0+=1
    return d0 
def nonRootWithOneChild(theList): # returns the number of non-root vertices with one child
    d1=0
    theList=childrenOfVertices(theList)
    theList.pop(0) #removes the first element of the list(children of the root)
    for item in theList:
        if item==1: #finds all the ones in the list
            d1+=1
    return d1
def returnVertices(myTree): #Write a function that returns the number of vertices in a plane tree.
    counter=1 #always adds itself
    for child in myTree:#goes through every child
        len(child)
        counter+=returnVertices(child)#how many vertices are at or below this vertex
    return counter
def returnEdges(myTree):
    n=returnVertices(myTree)-1
    return n
def energyFunction(mytree): #where E(T)= r + 2d0+ d1-n
    r= childrenOfRoot(mytree)
    d0= returnLeaves(mytree)
    d1= nonRootWithOneChild(mytree)
    n=returnEdges(mytree)

    return (-.4*r +(2.3*d0) +1.3*d1) - .1*n
#MARKOVCHAINMONTECARLO
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
def makeMove(i,j,dyckWord): #this function takes in a switch and outputs a new valid dyckWord
    counter1=0 # counts the ith positon of 1
    counter0=0 # counts the jth postion of 0
    oldWord=dyckWord
    newList=[] #list that contains dyckWords
    for letter in oldWord:
        if letter==1: 
            counter1+=1 #counting the 1s
            if counter1==i:
                newList.append(0)#change the ith position to 0
            else:
                newList.append(letter)
        if letter==0:
            counter0+=1 #counting the 0s
            if counter0==j: 
                newList.append(1) #change the jth position to 1
            else:
                newList.append(letter) #copy dyckWord into list
    if isValid(newList)== True: #isValid function checks if my 10 move makes a valid dyckWord
        return DyckWord(newList)
    else:
        return oldWord #if not a valid move
    
def markovMove(dyckWord):#markovMove function
    length=len(dyckWord) #pick random i,j between 1 and length
    i=randrange(1,length/2)
    j=randrange(1,length/2)
    newWord=makeMove(i,j,dyckWord)
    return newWord
def movingWithProb(currWord,newWord): #assigns probabilities to moves
    currWordEnergy=fasterEnergyFunction(currWord) # calculates energy
    newWordEnergy=fasterEnergyFunction(newWord) # then calculates energy
    probability= (min(1, exp(-newWordEnergy)/exp(-currWordEnergy))).n() #MCMC Part!!!!
    ran=random.random()
    if ran <= probability:  #sets new dyckword to be current state
        return newWord
    else: #keeps current dyckword as our current state
        return currWord
def myProject(startWord, mixingTimeT, sampleInterval, numOfSamples):# initialword, t , collect every x-amount of steps, number of samples that I want
    samples=[]#an empty list that will append the samples
    newWord=markovMove(startWord)
    currWord=startWord
    for i in range(mixingTimeT):  #I need my movingWithProb to run mixingTimeT amount of times (while loop)
        newWord=markovMove(currWord)
        currWord=movingWithProb(currWord,newWord)
    sampCount=0
    stepCount=0
    while sampCount < numOfSamples:  #I need my program to stop after I have collected numOfSamples amount of samples
        if stepCount==sampleInterval: #after sampleInterval amount of steps, append the currWord to my list 'samples' (for loop)
            samples.append(currWord)
            stepCount=0
            sampCount+=1
        else:
            newWord=markovMove(currWord)
            currWord=movingWithProb(currWord,newWord)
            stepCount+=1
    return samples #return samples
    


# In[1]:


# energy function rewritten by Anna to be faster and avoid recursion depth issues
def  fasterEnergyFunction(word): 
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
        #print(letter)
        #print(candidate_int_node_depths)
        #print(int_nodes)
         
    
    #print('root degree: ' + str(root_deg))
    #print('leaves: ' + str(num_leaves))
    #print('internal nodes: ' + str(int_nodes))
    #print('edges:  ' + str(num_edges))
        
        
    #r= childrenOfRoot(mytree)
    #d0= returnLeaves(mytree)
    #d1= nonRootWithOneChild(mytree)
    #n=returnEdges(mytree)

    return (-.4*root_deg +(2.3*num_leaves) +1.3*int_nodes) - .1*num_edges


# In[27]:


fasterEnergyFunction(DyckWord([1,1,1,0,1,0,0,1,1,1,0,1,0,0,0,1,0,0,1,1,1,0,0,0]))


# In[31]:


#When you want to run the program
outPutSamples= myProject(CompleteDyckWords_size(7).random_element(), 1000, 1000, 6)
for word in outPutSamples:
    tree=word.to_ordered_tree()
    tree.to_poset().show()


# In[14]:


import sys
import time
import timeit
start_time = time.time()

n=1000
startWord = DyckWord([1 for i in range(n)] + [0 for i in range(n)])
outPutSamples= myProject(startWord, 50000, 2000, 100)
#for word in outPutSamples:
#    tree=word.to_ordered_tree()
#    tree.to_poset().show()# ... do stuff


end_time = time.time()
print("Elapsed time was %g seconds" % (end_time - start_time))


#compute contact distances
from operator import add
cd_sums = [0 for i in range(2*n)]
for sample in outPutSamples:
    cds = contactDistances(sample)
    cd_sums = map(add, cd_sums, cds)

out1 = open('cds_ MCMC_thermo_1000n_50000ini_2000int.txt', 'w')
for d, s in enumerate(cd_sums):
    out1.write(str(d) + "\t" + str(s) + "\n")
out1.close()



out_file = open('MCMC_thermo_1000n_50000ini_2000int.txt', 'w')
for s in outPutSamples:
    out_file.write(str(s)+"\n")
    
out_file.close()


# In[17]:


list_plot(cd_sums[4:])


# In[30]:


def S(n,d): # cd for uniform distribution
    if d%2==1:
        return 0
    return 1 / (d/2 +1) * binomial(d, d/2) * binomial(2*n - d -1, n - d/2 -1)

p1 =list_plot(cd_sums[4:], color='red', size=5)
r = 100 / catalan_number(1000)
p2 = list_plot([r*S(1000, d) for d in range(4, 2000, 1)], color='green', size=5)
show(p1+p2)
print [r*S(1000, d).n() for d in range(0, 20, 1)]


# In[11]:


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


# In[12]:


test_word = DyckWord([1,1,1,0,1,0,0,1,0,0])
contactDistances(test_word)


# In[ ]:




