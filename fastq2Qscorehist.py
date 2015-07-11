#!/usr/bin/python3
#
# fastq2Qscorehist.py
# 
import matplotlib.pyplot as plt

"""read a fastq file and plot histogram of quality scores"""

def readFastq(filename):
    """read fastq file and returns sequences and qualityscores"""
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities
    
def phred33ToQ(qual):
    "converts phred33 values to quality score"""
    return ord(qual) - 33

def createHist(qualityStrings):
    """helper function to create histogram of quality scores"""
    hist = [0]*50
    for read in qualityStrings:
        for phred in read:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist

if __name__=="__main__":
    
    #read test fastq file and create histogram of quality scores
    seq, qual = readFastq("testfile.fastq")   
    h = createHist(qual)
    plt.bar(range(len(h)), h)
    plt.show()
