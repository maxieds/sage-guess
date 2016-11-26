#### GuessConfig.py : 
#### Constants and other definitions for the guess package bundle
#### Author:  Maxie D. Schmidt
#### Created: 2016.11.25

from ColorizeString import print_color_string

DEBUGGING = True
DEFAULT_NUM_PRIMES = 64
DEFAULT_NUM_ROWS = 32
MAX_OEIS_RESULTS = 2
FIRST, LAST, NULL = 0, 1, None
X, Y = 0, 1

def TODO(value): return value

def print_debug(msg, debugging = DEBUGGING, fg = 'FGWHITE', bg = 'BGRED'): 
     if not debugging: 
          return None
     print_color_string("!!! DEBUG !!!", opts = [fg, bg, 'BOLD', 'UL'])
     print_color_string(" : ", opts = ['NORMAL'])
     print msg
## 

def IsZero(elem): 
     if isinstance(elem, int) or isinstance(elem, long): 
          return elem == 0
     else: 
          return elem.is_zero()
     ##
##

def IsOne(elem): 
     if isinstance(elem, int) or isinstance(elem, long): 
          return elem == 1
     else: 
          return elem.is_one()
     ##
##

