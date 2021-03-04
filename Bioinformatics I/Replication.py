# INPUT: Text (Dna string), Pattern (shorter Dna fragment being compared to Text)
# slides window down the Dna string and compares it to a given pattern
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i : i + len(Pattern)] == Pattern:
            count = count + 1
    return count

# INPUT: Text (Dna string), k (length of kmer, or length of patterns we're looking for)
def FrequentWords(Text, k):
    words = [] # list of most common appearing kmers
    freq = FrequencyMap(Text, k) 
    # freq is a dictionary whose keys are kmers and values are the number of times the kmer appears 
    m = max(freq.values()) # finds key with highest value, AKA most frequent kmer

    # appends all kmers that appear m times into the words[] list
    for key in freq: 
        if (freq[key] == m):
            words.append(key)
    return words

def FrequencyMap(Text, k):
    freq = {}
    # overall runtime of FrequencyMap is O(n^2)

    for i in range(len(Text) - k + 1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0 # sets each kmer to a key whose default value is 0

        # iterates through whole DNA string again, comparing segments to the kmer key
        for i in range(len(Text) - k + 1):
            if (Pattern == Text[i : i + k]):
                freq[Pattern] = freq[Pattern] + 1
    return freq

def ReverseComplement(Pattern):
    reverseComplement = ""
    for i in range(len(Pattern)):
        if (Pattern[i] == 'A'):
            reverseComplement += 'T'
        elif (Pattern[i] == 'T'):
            reverseComplement += 'A'
        elif (Pattern[i] == 'G'):
            reverseComplement += 'C'
        else:
            reverseComplement += 'G'
    
    return reverseComplement[::-1] # reverses string 

# INPUT: two strings, the genome (named Text) and Pattern
# OUPUT: integer list specifying starting positions where Pattern is a substring of Genome
def PatternMatching(Text, Pattern):
    positions = []
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i : i + len(Pattern)] == Pattern:
            positions.append(i)
    return positions



Genome = ""
Pattern = "CTTGATCAT "

data = PatternMatching(Genome, Pattern)
print(*data, sep='  ') # use this snippet here to omit commas when printing to console