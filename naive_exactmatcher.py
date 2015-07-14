#!/usr/bin/python3
#
# naive_exactmatcher.py
#
# Usage: naive_exactmatcher.py <scafoldfilename> <fastqfilename>
#
"""A script implementing the naive exact pattern matching algorithm to match 
sequenced reads to a scafold

File containg single scaffold and fastq file containg reads must be supplied 
at command line. Shorter readlength can be evaluated by providing the value of 
inp. Deafult is full read length.
"""

inp = None #set value if evaluating shorter input(eg, 30)


import sys

def readscaffold(filename):
    """parse a file containing a scaffold and return one string of nuc sequence"""
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome



def reverseComplement(s):
    """return a reverse complement of a DNA nucleotide string"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    newseq = ''
    for base in s:
        newseq = complement[base] + newseq
    return newseq

def readFastq(filename):
    """parse a fastq file and return a list of sequences and quality scores"""
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def naive_exactmacth(readseq, genomeseq):
    """function implementing the naive pattern matching algorithm"""
    occurrences = []
    for i in range(len(genomeseq) - len(readseq) + 1):
        match = True
        for j in range(len(readseq)):
            if genomeseq[i+j] != readseq[j]:
                match = False
                break
        if match:
          occurrences.append(i)
    return occurrences
    
if __name__=="__main__":
    
    
    try:
        scaf = sys.argv[1]
        readfile = sys.argv[2]
        scaffold = readscaffold(scaf)
        reads = readFastq(readfile)[0]
        numMatched = 0
        n = 0
        
       #loop through reads and identify exact matches to scaffold
        for r in reads:
              r = r[:inp]
              matches = naive_exactmacth(r, scaffold)
        
          # return only one match if read and reverse complements are palindromes
              if reverseComplement(r) != r: 
                matches.extend(naive_exactmacth(reverseComplement(r), scaffold))        
              n += 1
              if len(matches) > 0:
                  numMatched += 1
        print('%d / %d reads matched the genome exactly!' % (numMatched, n))
    except :
             print("Wrong usage type: (naive_exactmatcher.py scafoldfilename fastqfilename)")    
