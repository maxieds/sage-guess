
r"""
GuessConfig.py 

Constants and other definitions for the guess package bundle

AUTHORS: 
- Maxie D. Schmidt (Created: 2016.11.25) 

"""

from ColorizeString import print_color_string

DEBUGGING = True
DEFAULT_NUM_PRIMES = 64
DEFAULT_NUM_ROWS = 32
MAX_OEIS_RESULTS = 2
FIRST, LAST, NULL = 0, 1, None
X, Y = 0, 1

def TODO(value): 
     r"""
     Returns its input. Primarily used for documenting features / code that 
     needs to be updated in the source. 

     INPUTS: 
     - ``value`` -- Any value. 
     """
     return value
## 

def IsZero(elem): 
     r""" 
     Returns whether the input element represents zero. 

     INPUTS: 
     - ``elem`` -- An integer-like Python or Sage variable
     """ 
     if isinstance(elem, int) or isinstance(elem, long): 
          return elem == 0
     else: 
          return elem.is_zero()
     ##
##

def IsOne(elem): 
     r""" 
     Returns whether the input element represents one. 

     INPUTS: 
     - ``elem`` -- An integer-like Python or Sage variable
     """ 
     if isinstance(elem, int) or isinstance(elem, long): 
          return elem == 1
     else: 
          return elem.is_one()
     ##
##

def IsInteger(elem): 
     r""" 
     Returns whether the input element is an integer of some type. 

     INPUTS: 
     - ``elem`` -- An integer-like Python or Sage variable
     """ 
     if isinstance(elem, int) or isinstance(elem, long): 
          return True
     else: 
          return elem.is_integer()
     ##
##

