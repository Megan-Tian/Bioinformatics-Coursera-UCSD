# Which DNA Patterns Play The Role of Molecular Clocks? (Part 1 + 2)
# Week 3 and 4 
# Biology Meets Programming: Bioinformatics for Beginners by University of California San Diego

import random 
import MotifAlgorithms

###----------SEARCH FUNCTIONS (Gibbs, Random, Greedy)----------###

# GibbsSampler function 
def GibbsSampler(Dna, k, t, N):
    BestMotifs = [] # output variable
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(1,N):
        i = random.randint(0,t-1)
        ReducedMotifs = []

        for j in range(0,t):
            if j != i:
                ReducedMotifs.append(Motifs[j])

        Profile = ProfileWithPseudocounts(ReducedMotifs)
        Motifs[i] = ProfileGeneratedString(Dna[i], Profile, k)

        if (Score(Motifs) < Score(BestMotifs)):
                BestMotifs = Motifs

    return BestMotifs 

# Detailed explanation of GreedyMotifSearch: https://www.mrgraeme.com/greedy-motif-search/
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    
    for i in range(len(Dna[0])-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])

        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
            
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
        
    return BestMotifs

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = [] # output variable
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 

###----------SUBROUTINES FOR SEARCH FUNCTIONS----------###
# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    normalProb = {}
    probSum = 0
    for key in Probabilities:
        probSum += Probabilities[key]
    
    for key in Probabilities:
        normalProb[key] = Probabilities[key]/probSum
        
    return normalProb

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
        
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    kmer = '' # output variable
    p = random.uniform(0, 1)
    for i in Probabilities:
        p -= Probabilities[i]
        if (p <= 0):
            kmer = i
            return kmer

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile) - Pr = probability func
def Pr(Text, Profile):
    # insert your code here
    p = 1 # p = probability variable
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    
    return p

# Input: A string Text, an integer k, and a 4 x k matrix Profile.
# Output: A Profile-most probable k-mer in Text.
def ProfileMostProbableKmer(Text, k, Profile):
    Kmer = ""
    maxProb = -1
    for i in range(len(Text)-k+1):
        tempProb = Pr(Text[i:i+k], Profile) 
        if tempProb > maxProb:
            maxProb = tempProb
            Kmer = Text[i:i+k]
    return Kmer  

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers. We can compute Score(Motifs) by first constructing Consensus(Motifs) and then summing the number of symbols in the j-th column of Motifs that do not match the symbol in position j of the consensus string
def Score(Motifs):
    consensus = Consensus(Motifs)
    count = 0
    for motif in Motifs:
        for index in range(len(motif)):
            if motif[index] != consensus[index]:
                count += 1
    return count

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)    
    consensus = ""
    
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
   
    return consensus 

# Input: list of kmer Motifs
# Output: profile matrix of motifs as a dict of lists. Build a Profile matrix by dividing each element of the count matrix by the amount of rows in the Motifs matrix

def Profile(Motifs):
    t = len(Motifs) # number of rows
    k = len(Motifs[0]) # number of columns
    profile = Count(Motifs)
    
    for i in profile:
        for j in range(k):
            profile[i][j] = profile[i][j] / t
    
    return profile

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs) # number of rows
    k = len(Motifs[0]) # number of columns
    profile = CountWithPseudocounts(Motifs)
    
    for i in profile:
        for j in range(k):
            profile[i][j] = profile[i][j] / (t + 4)
    
    return profile

# Input:  Motifs, a list of kmers (strings)
# Output: Count(Motifs), the count matrix of Motifs as a dictionary of lists
def Count(Motifs):
    count = {} # initializing the count matrix
    
    k = len(Motifs[0]) # length of the first kmer, number of columns in Motifs array, length of each list in the dictionary (all kmers are same length)
    
    for symbol in "ACGT": # symbol is arbitrary counter like "i"
        count[symbol] = [] # count matrix now has keys A, C, T, and G all with values of empty list
        for j in range(k):
            count[symbol].append(0) # count matrix now has keys A, C, G, and T all with values of a list of zeroes of length equal to the length of a kmer
            
    t = len(Motifs) # length of Motifs, a list of kmers (strings)
    
    for i in range(t): # for each kmer in Motifs
        for j in range(k): # for each element of the kmer
            symbol = Motifs[i][j] # assigns the key (symbol) to a nucleotide (ACGT) in Motifs
            
            #count[symbol] corresponds to the key of the dictionary count
            #count[symbol][j] corresponds to the position in the list assigned to the key
            count[symbol][j] += 1 # adds 1 to the position in the list assigned to the key

    return count

def CountWithPseudocounts(Motifs):
    count = {} 
    k = len(Motifs[0]) 
    
    for symbol in "ACGT": 
        count[symbol] = [] 
        for j in range(k):
            count[symbol].append(1) # only change from Count(Motifs) is to set the default initial matrix values to 1 instead of 0
            
    t = len(Motifs) 
    
    for i in range(t): 
        for j in range(k): 
            symbol = Motifs[i][j]
            count[symbol][j] += 1 
    
    return count

# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t). Choose a random k-mer from each of t different strings Dna, and returns a list of t strings.
def RandomMotifs(Dna, k, t):
    randomMotifs = []
    for i in range(len(Dna)):
        r = random.randint(1, len(Dna[0]) - k)
        randomMotifs.append(Dna[i][r : r + k])
    return randomMotifs

def Motifs(Profile, Dna):
    motifs = []
    for i in range(len(Dna)):
        m = ProfileMostProbableKmer(Dna[i], len(Profile), Profile)
        motifs.append(m)
    return motifs

