# set Text equal to the Vibrio cholerae oriC and k equal to 10
oriC = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
k = 10

# print the result of calling FrequentWords on Text and k. Returns most frequent strings of length k from the Dna string
print (FrequentWords(oriC, k))

# TODO - figure out how to read oriC_Vibrio_Cholerae text file - course does auto-import with built in code editor which can't be replicated here

#----------------------------------------------------------#

# Input:  Strings Pattern and Text along with an integer d
# Pattern = the longer string we're finding d-mismatched k-mers in 
# Text = the shorter string we're comparing to the longer one in iterations
# Output: A list containing all starting positions where Pattern appears as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i : i + len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q, or number of mismatches between p and q
def HammingDistance(p, q):
    # assuming p and q are the same length
    hDistance = 0
    for i in range (len(p)):
        if (p[i] != q[i]):
            hDistance += 1
    return hDistance

# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
    positions = [] # output variable
    s = SkewArray(Genome)
    minSkew = min(s)
    for i in range (len(s)):
        if (s[i] == minSkew):
            positions.append(i)
            # adding the position of the skew not the value, so append(i) not append(s[i])
    return positions

def SkewArray(Genome):
    skew = [0]
    score = {"A":0, "T":0, "C":-1, "G":1}
    for i in range(1,len(Genome)+1):
            skew.append(score[Genome[i-1]] + skew[i-1])
    return skew

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    for i in range (0, len(Genome)):
        if (Genome[i:i+len(Pattern)] == Pattern):
            positions.append(i)
    return positions

# Input:  A DNA string Pattern
# Output: The reverse complement of Pattern
def ReverseComplement(Pattern):   
    return Reverse(Complement(Pattern))

def Complement(Pattern):
    complement = ""
    for i in range(0, len(Pattern)):
        if (Pattern[i] == 'A'):
            complement += 'T'
        elif (Pattern[i] == 'T'):
            complement += 'A'
        elif (Pattern[i] == 'G'):
            complement += 'C'
        elif (Pattern[i] == 'C'):
            complement += 'G'
    return complement 

# Input:  A string Pattern
# Output: The reverse of Pattern
def Reverse(Pattern):
    reverse = Pattern[::-1]
    return reverse


# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if (freq[key] == m):
            words.append(key)
    return words

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    # hint: your code goes here!    
        for i in range(n-k+1):
            if (Pattern == Text[i:i+k]):
                freq[Pattern] = freq[Pattern] + 1
    return freq

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count



