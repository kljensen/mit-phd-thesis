#!/usr/bin/env python
import string
import sys
import math

# Usage
def usage():
    print sys.argv[0], ": make a weight matrix from an alignment file"
    print "Usage: ", sys.argv[0], "<alignment file>"



# a PWM class, basically just a list of dictionaries
# with a few methods for building/accessing
class Pwm:
    """A position weight matrix class"""
    def __init__(self, aln):

        #  find out what characters occurs in seqs
        self.chars = self.getChars(aln)
        self.length = aln.support()
        self.width = aln.wid()

        # initialize list of dicts
        # cm = count matrix
        # fm = frequency matrix
        self.cm = []
        self.fm = []
        for j in range(aln.width):
            self.cm.append( dict([(char, 0)for char in self.chars]) )
            self.fm.append( dict([(char, 0.0)for char in self.chars]) )

        # fill in the counts
        for j in range(self.wid()):
            for i in range(self.len()):
                c = aln.getChar(i,j)
                self.cm[j][c] = self.cm[j][c] + 1

        # calculate the frequency matrix
        for c in self.chars:
            for j in range(self.wid()):
                self.fm[j][c] = float(self.cm[j][c]) / float(self.len())

        self.calcEntropy()

    def len(self):
        return self.length
    def wid(self):
        return self.width
    def getChars(self, aln):
        chars = []
        for j in range(aln.width):
            for i in range(aln.support()):
                c = aln.getChar(i,j)
                if chars.count(c) == 0:
                    chars.append(aln.getChar(i,j))
        return chars
    def getFreq(self, pos, char):
        return self.fm[pos][char]
    def getCount(self, pos, char):
        return self.cm[pos][char]
    def printPwm(self):
        print "Frequency Matrix:"
        for char in self.chars:
            s = "%s" % (char)
            for pos in range(self.wid()):
                s = s + (" %1.3f" % ( self.getFreq(pos, char) ))
            print s

    def printPwmDNA(self):
        print "Frequency Matrix:"
        for char in ['A', 'T', 'G', 'C']:
            s = "%s" % (char)
            for pos in range(self.wid()):
                s = s + (" %1.3f" % ( self.getFreq(pos, char) ))
            print s
        print ""

    # returns a vector, each member of which is
    # the entropy at a pos in the pwm
    def calcEntropy(self):
        self.entropy = []
        base = 2
        for pos in range(self.wid()):
            self.entropy.append(0)
            for char in self.chars:
                freq = self.getFreq(pos, char)
                if (freq > 0):
                    self.entropy[pos] = self.entropy[pos] \
                            - freq * math.log(freq, base)

    def getEntropy(self, pos):
        return self.entropy[pos]

    # returns a vector, each member of which is
    # the info content at a pos in the pwm
    # take a dictionary containing background 
    # frequencies of various characters
    def calcInfo(self, bg):
        bgInfo = 0.0
        base = 2
        for char, prior in bg.iteritems():
            if prior > 0:
                bgInfo = bgInfo - prior * math.log(prior,base)

        self.info = []
        self.ic = 0
        for pos in range(self.wid()):
            self.info.append(bgInfo - self.getEntropy(pos))
            self.ic = self.ic + bgInfo - self.getEntropy(pos)


    def getInfo(self, pos):
        return self.info[pos]
    def getIC(self):
        return self.ic

    def printInfo(self):
        for pos in range(self.wid()):
            print "Position %d: entropy = %1.3f, information = %1.3f" \
                    % (pos+1, self.getEntropy(pos), self.getInfo(pos))
        print "\nTotal Information Content = %1.3f" % (self.getIC())

    # compute a bitScore match of the pwm against a sequence
    def score(self, seq, bg):
        total = 0
        for pos in range(self.wid()):
            char = seq[pos]
            #  potential need for error checking here!
            freq = self.getFreq(pos, char)
            prior = bg[char]
            if freq != 0 and prior != 0:
                total = total + math.log(freq/prior,2)
        return total





# A simple sequence alignment class
class Alignment:
    """A simple sequence alignment class"""
    def __init__(self, label):
        self.label = label
        self.seqs = []
        self.width = 0
    def addSeq(self, seq):
        self.seqs.append(seq)
        if(len(seq) > self.width):
            self.width = len(seq)
    def printAlignment(self):
        for seq in self.seqs:
            print seq
    def support(self):
        return len(self.seqs)
    def wid(self):
        return self.width
    def getChar(self, seq, pos):
        return self.seqs[seq][pos]
    def makePWM(self):
        self.pwm = Pwm(self)
    def printPWM(self):
        self.pwm.printPwm()
    def printPWMDNA(self):
        self.pwm.printPwmDNA()
    def calcInfo(self, bg):
        self.pwm.calcInfo(bg)
    def printInfo(self):
        self.pwm.printInfo()
    def selfScore(self, bg):
        mean = 0
        for seq in self.seqs:
            score = self.pwm.score(seq, bg)
            print "Sequence %s has score s = %0.3f"  % (seq, score)
            mean = mean + score
        mean = mean / self.support()
        print "Mean score = %.3f" % (mean)



def readAlignFile(alignFH, label):
    aln = Alignment(label)
    for line in alignFH:
        line = string.strip(line)
        aln.addSeq(line)
    return aln

# Main
def main():

    # see if we got the right number of command line args
    if len(sys.argv) != 2:
        usage()
        sys.exit(2)

    # get command line paramters
    fname1 = sys.argv[1]

    # open the alignment file
    try:
        alignFH = open(fname1, "r")
    except IOError, (errno, strerror):
        print "Error %s: %s" % (errno, strerror)
        sys.exit()


    aln = readAlignFile(alignFH, "My alignment")

    # close alignment file
    try:
        alignFH.close()
    except IOError, (errno, strerror):
        print "Error %s: %s" % (errno, strerror)
        sys.exit()

    # do all our fancy shizzle
    #aln.printAlignment()
    aln.makePWM()
    aln.printPWMDNA()
    aln.calcInfo({'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25})
    aln.printInfo()
    print ""
    aln.selfScore({'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25})

# execute main
if __name__ == "__main__":
    sys.exit(main())


