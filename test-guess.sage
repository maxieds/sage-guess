attach("ColorizeString.py")
attach("GuessingFormulasForSequences.py")
attach("Guess.py")
attach("SearchSequenceFactors.py")
attach("SequenceElement.py")
attach("SpecialSequences.py")

from Guess import *
from SpecialSequences import *
from ColorizeString import print_color_string
from SearchSequenceFactors import *

def print_match_summary(fmatches): 
     for (fidx, fmatch) in enumerate(fmatches): 
          header_string = " ==== Formula Match # %d / %d ==== " % \
               (fidx + 1, len(fmatches))
          print_color_string(header_string, 
                             opts = ['FGWHITE', 'BGBLACK', 'BOLD', 'UL'])
          print ""
          fmatch.print_summary()
     ##
     print "\n"
## 

#sm = guess([0, 1, 4, 9, 16, 25, 36])
#print_match_summary(sm)

#sm = guess([factorial(n+1) * harmonic_number(n+1, 1) for n in range(0, 8)])
#print_match_summary(sm)

#seq = [factorial(2*n+1) * harmonic_number(2*n+1, 1) for n in range(0, 8)]
#sm = guess(seq)
#print_match_summary(sm)

#seq = [S1(3*n+1, 2) for n in range(0, 6)]
#print seq
#sm = guess(seq)
#print_match_summary(sm)

#seq = [S1(3*n+1, 2*n+2) for n in range(1, 8)]
#print seq
#sm = guess(seq)
#print_match_summary(sm)

seq = [13 * S1(3*n+1, 2*n+2) for n in range(1, 8)]
print seq
sm = guess(seq)
print_match_summary(sm)

seq = [13 * n * S1(3*n+1, 2*n+2) for n in range(1, 8)]
print seq
sm = guess(seq)
print_match_summary(sm)


