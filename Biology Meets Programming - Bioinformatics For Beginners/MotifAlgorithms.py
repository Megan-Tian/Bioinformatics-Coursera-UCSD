import Motifs

###----------DATASET AND OTHER PARAMETERS----------###
# 10 strings in DosR (tuberculosis) dataset
Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",      
"CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",        
"ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",        
"GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",        
"GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",        
"CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",        
"GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",        
"GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",        
"GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",       
"TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

k = 15 # length of the kmer we're looking for
t = 10 # number of strings in Dna
N = 100 # iterations in a cycle of ex. GibbsSampler, Randomized Motif, etc.


###----------GIBBS SAMPLER ON THE DATASET----------###
# Call GibbsSampler(Dna, k, t, N) 20 times and store the best output in a variable
BestGibbsMotifs = GibbsSampler(Dna, k, t, N)

for i in range(20):
    m = GibbsSampler(Dna, k, t, N)
    if Score(m) < Score(BestGibbsMotifs):
        BestGibbsMotifs = m

print("GIBBS SAMPLER SEARCH")
print("Kmer length = " + k)
print("The lowest scoring motif(s): ")
print(BestGibbsMotifs)
print("The score of these motifs is: ")
print(Score(BestGibbsMotifs))


###----------GREEDY MOTIF SEARCH ON THE DATASET----------###
BestGreedyMotifs = GreedyMotifSearch(Dna, k, t)

print("GREEDY SEARCH")
print("Kmer length = " + k)
print("The lowest scoring motif(s): ")
print(BestGreedyMotifs)
print("The score of these motifs is: ")
print(Score(BestGreedyMotifs))


###----------GREEDY MOTIF SEARCH WITH PSEUDOCOUNTS ON THE DATASET----------###
BestGibbsWithPseudocountsMotifs = GreedyMotifSearchWithPseudocounts(Dna, k, t)

print("GREEDY WITH PSEUDOCOUNTS SEARCH")
print("Kmer length = " + k)
print("The lowest scoring motif(s): ")
print(BestGibbsWithPseudocountsMotifs)
print("The score of these motifs is: ")
print(Score(BestGibbsWithPseudocountsMotifs))

###----------RANDOMIZED MOTIF SEARCH ON THE DATASET----------###
# Call RandomizedMotifSearch(Dna, k, t) N times, storing the best-scoring set of motifs resulting from this algorithm in a variable 
BestRandomizedMotifs = RandomizedMotifSearch(Dna, k, t)

for i in range(N):
    m = RandomizedMotifSearch(Dna, k, t)
    if Score(m) < Score(BestRandomizedMotifs):
        BestRandomizedMotifs = m

print("GIBBS SAMPLER SEARCH")
print("Kmer length = " + k)
print("The lowest scoring motif(s): ")
print(BestRandomizedMotifs)
print(Score(BestRandomizedMotifs))

###----------RANDOMIZED MOTIF SEARCH WITH PSEUDOCOUNTS ON THE DATASET (not included in the course) ----------###
