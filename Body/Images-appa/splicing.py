#!/usr/bin/env python
import string
import sys
import re

# Usage
def usage():
    print sys.argv[0], ": searches for IUPAC regular expressions in a fasta file"
    print "Usage: ", sys.argv[0], "<motif file> <fasta file>"

class IntronSeq:
    """A simple intron class"""
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop
        self.length = stop-start
    def Start(self):
        return self.start
    def Stop(self):
        return self.stop
    def Length(self):
        return self.length

# A simple fasta sequence class
class FastaSeq:
    """A simple fasta sequence class"""
    def __init__(self, label, seq):
        self.label = label
        self.seq = seq
    def getSeq(self):
        return self.seq
    def getLabel(self):
        return self.label
    def appendLabel(self, newLabel):
        self.label = self.label + newLabel
    def appendSeq(self, newSeq):
        self.seq = self.seq + newSeq
    def getLength(self):
        return len(self.seq)
    def getSubSeq(self, start, length):
        return self.seq[start:start+length]
    def getSubSeq2(self, start, stop):
        return self.seq[start:stop]

    # finds a regular expression in the sequence and returns
    # a list of integers specifying where the pattern starts
    def findRegEx(self, pattern):
        lastOffset = 0
        count = 0
        myRegEx =  re.compile(pattern)
        positions = []
        while 1:
            match = myRegEx.search(self.seq, lastOffset)
            if match:
                count = count + 1
                lastOffset = match.start() + 1
                positions.append( match.start() )
            else:
                break
        return positions

    # finds a regular expression in the sequence and returns
    # a list of integers specifying where the pattern ends
    def findRegExEnd(self, pattern):
        lastOffset = 0
        count = 0
        myRegEx =  re.compile(pattern)
        positions = []
        while 1:
            match = myRegEx.search(self.seq, lastOffset)
            if match:
                count = count + 1
                lastOffset = match.end() + 1
                positions.append( match.end() )
            else:
                break
        return positions

def getIntrons(mySeq, fpRegEx, tpRegEx):
    fpSites = mySeq.findRegEx(fpRegEx)
    tpSites = mySeq.findRegExEnd(tpRegEx)
    i=0
    j=0
    introns = []
    while i<len(fpSites):
        while j<len(tpSites) and tpSites[j] <= fpSites[i]:
            j=j+1
        #print "Intron from %d to %d!" % (fpSites[i], tpSites[j])
        newIntron = IntronSeq(fpSites[i], tpSites[j])
        introns.append(newIntron)
        while i<len(fpSites) and fpSites[i] <= tpSites[j]:
            i=i+1
    return introns


# read in a fasta file and return a list of fastaSeq objects
def readFastaFile(fh):
    numSeq = 0
    mySeqs = []

    #  read in a line
    for line in fh:

        # strip off white space
        line = string.strip(line)

        if line[0] == ">":
            # we hit a label line
            newSeq = FastaSeq(line[1:], "")
            mySeqs.append( newSeq )
            numSeq = numSeq + 1
        else:
            # not a label line
            mySeqs[ numSeq - 1 ].appendSeq(line)

    return mySeqs

def iupac2regex(motif):
    regex = motif
    i2rTable = {    'N': '[ATGC]',
                    'K': '[GT]',
                    'B': '[CGT]',
                    'S': '[CG]',
                    'V': '[ACG]',
                    'Y': '[CT]',
                    'W': '[AT]',
                    'D': '[AGT]',
                    'R': '[GA]',
                    'M': '[AC]',
                    'H': '[ACT]'}
    for old,new in i2rTable.iteritems():
        regex = string.replace(regex, old, new)
    return regex

# Main
def main():

    # see if we got the right number of command line args
    if len(sys.argv) != 3:
        usage()
        sys.exit(2)

    # get command line paramters
    fname1 = sys.argv[1]
    fname2 = sys.argv[2]

    # open the motif file
    try:
        motifFH = open(fname1, "r")
    except IOError, (errno, strerror):
        print "Error %s: %s" % (errno, strerror)
        sys.exit()

    # open the fasta file
    try:
        fastaFH = open(fname2, "r")
    except IOError, (errno, strerror):
        print "Error %s: %s" % (errno, strerror)
        sys.exit()

    # read in the DNA sequences
    myDnaSeqs = readFastaFile(fastaFH)

    # read in our regexs
    regexs = []
    patterns = []
    for line in motifFH:
        regexs.append( iupac2regex( string.strip(line) ) )
        patterns.append(string.strip(line))
    # assume 0 = 5' && 1 = 3'

    # splice each seq
    for seq in myDnaSeqs:
        introns = getIntrons(seq, regexs[0], regexs[1])
        splicedSeq = ""
        last = 0
        for intron in introns:
            #if intron.Length() >= 300 and intron.Length() <= 600:
            print "found 5' at %d" % (intron.Start())
            print "found 3' at %d" % (intron.Stop()-len(patterns[1]))
            splicedSeq = splicedSeq + seq.seq[last:intron.Start()]
            last = intron.Stop()
        splicedSeq = splicedSeq + seq.seq[last:seq.getLength()]
        print "Concatination of exons:"
        print splicedSeq


# execute main
if __name__ == "__main__":
    sys.exit(main())


