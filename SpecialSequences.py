
r""" 
SpecialSequences.py 

Implements several well-known triangular sequences in two indices 
and a couple of other special function forms still not stock in Sage. 

AUTHORS: 
- Maxie D. Schmidt (Created: 2016.05.19) 

"""

import sys
import numpy as np
from sage.all import *
from sage.symbolic.function_factory import function_factory
from SequenceElement import SequenceGenerator
from GuessConfig import *

def pochhammer(x, n):
     r""" 
     Implements the (ordinary) Pochhammer symbol, 
     $(x)_n = x(x+1)(x+2) ... (x+n-1)$. 

     INPUTS: 
     - ``x`` -- The indeterminate parameter x
     - ``n`` -- The function index 

     EXAMPLES: 

     We can use the Pochhammer symbol to implement the double factorial function, 
     (2n-1)!! = 2^{n} (1/2)_n, for n >= 1: 

     sage: from SpecialSequences import pochhammer
     sage: fact2 = lambda n: int((2 ** n) * pochhammer(0.5, n))
     sage: fact2_seq = [fact2(2) for n in range(1, 8)]
     sage: print fact2_seq
     [1, 3, 15, 105, 945, 10395, 135135]
     sage: from Guess import *
     sage: seqnums = guess_oeis_formula(fact2_seq)
     sage: seqnums[0].print_summary()
        >> A001147: Double factorial of odd numbers: a(n) = (2*n-1)!! = 1*3*5*...*(2*n-1).
        >> First Sequence Values:  [1, 1, 3, 15, 105, 945, 10395, 135135, 2027025, 34459425]

     """
     product_factor_func = lambda i: x + i
     factor_list = map(product_factor_func, range(0, n))
     return prod(factor_list)
## 

def harmonic_number(n, r = 1):
     r""" 
     Implements the r-order harmonic numbers, $H_n^{(r)}$, which form the 
     sequence of partial sums of the Riemann zeta function, $\zeta(r)$, 
     for integers n >= 0. 

     INPUTS: 
     - ``n`` -- The sequence index 
     - ``r`` -- The order of the terms in the partial sums 
     """
     k = var('k')
     return sum(k ** -r, k, 1, n)
## 

def identity_func(): 
     r""" 
     Implements a one-dimensional identity function. 
     """
     return lambda n: n
##

def zero_func(): 
     r"""
     Implements a one-dimensional zero function. 
     """
     return lambda n: 0
## 

def one_func(): 
     r"""
     Implements a one and two-dimensional integer 1-valued unit function. 
     """ 
     return lambda n, k = 1: 1
## 

def constant_func(constant): 
     r""" 
     Implements a one-dimensional constant function. 

     INPUTS: 
     - ``constant`` -- A scalar-valued constant 
     """
     return lambda n: constant
## 

def sequence_func(seq_data):
     r"""
     Implenents a single-input function that returns values from a fixed 
     input sequence starting from the index n >= 0. 

     INPUTS: 
     - ``seq_data`` -- A finite list of sequence terms 
     """
     return lambda n: seq_data[n] if n >= 0 and n < len(seq_data) else 0

class TriangularSequence(object):
     r""" 
     TriangularSequence : 
     Implements the general form of a two-index, six-parameter 
     triangular combinatorial sequence defined recursivly in 
     Exercise #6.94 of Concrete Mathematics: 
          T(n, k) = (an+bk+g) T(n-1, k) + (a'n+b'k+g') T(n-1,k-1) + [n = k = 0]. 
     Special cases of this sequence form include the binomial coefficients 
     (in the simplest case), the (unsigned) Stirling numbers of the first and 
     second kinds, and the first-order and second-order Eulerian number triangles. 
     Other variants form the coefficients in the symbolic polynomial 
     expansions of generalized factorial functions of the form 
          R(R+a)(R+2a) ... (R+(n-1)a), and 
          b(a+b)(2a+b) ... (a(n-1)+b), 
     for indeterminates R, and scalar-valued a, b (not both zero). 
     """

     def __init__(self, a, b, g, ap, bp, gp, max_rows = DEFAULT_NUM_ROWS):
          r""" 
          Initializes the TriangularSequence object based on its 
          (integer-valued) parameters. 

          INPUTS: 
          - ``a`` -- Recurrence parameter (typically integer-valued) 
          - ``b`` -- Recurrence parameter (typically integer-valued)
          - ``g`` -- Recurrence parameter (typically integer-valued)
          - ``ap`` -- Recurrence parameter (typically integer-valued) 
          - ``bp`` -- Recurrence parameter (typically integer-valued) 
          - ``gp`` -- Recurrence parameter (typically integer-valued)
          - ``max_rows`` -- The maximum number of stored triangle rows to 
                            compute at initialization of the object 
                            (for quicker repeated lookups) 
          """
          self.alpha = a 
          self.beta = b
          self.gamma = g
          self.alpha_prime = ap
          self.beta_prime = bp
          self.gamma_prime = gp
          self.stored_rows = max_rows
          self.rec_data = []
          self.generate_rows(1, max_rows)
     ## 

     def generate_rows(self, row_min, row_max):
          r"""
          Generates the stored triangle data for rows indexed by the 
          upper input n = 0, 1, ...., row_max. 

          INPUTS: 
          - ``row_min`` -- Non-negative integer specifying the minimum row 
                           index to compute. Useful for rebuilding, or 
                           adding new rows to an existing table. 
          - ``row_max`` -- Non-negative integer specifying the number of 
                           rows to compute 
          """
          [a, b, g, ap, bp, gp] = [self.alpha, self.beta, self.gamma, \
                                   self.alpha_prime, self.beta_prime, \
                                   self.gamma_prime]
          for n in range(row_min, row_max + 1): 
               row_data = []
               for k in range(0, n + 1): 
                    recnk = (a * n + b * k + g) * self.get_data(n - 1, k) + \
                            (ap * n + bp * k + gp) * self.get_data(n - 1, k - 1)
                    row_data += [recnk]
               ##
               self.rec_data += [row_data]
          ##
          self.stored_rows = row_max
     ## 

     def get_data(self, n, k): 
          r"""
          Returns the triangle element at row n and column k. 

          INPUTS: 
          - ``n`` -- The row index of the element to be returned 
          - ``k`` -- The lower column index of the element to be returned 
          """
          if n < 0 or k < 0 or k > n or n > self.stored_rows:
               return 0
          elif n == 0 and k == 0: 
               return 1
          elif n > self.stored_rows:
               self.generate_rows(self.stored_rows + 1, n)
          ##
          return self.rec_data[n][k]
     ## 

     def print_table_rows(self, row_min, row_max):
          r"""
          Prints the triangle data stored locally from row row_min to 
          row index row_max. 

          INPUTS: 
          - ``row_min`` -- The lower row index to print
          - ``row_max`` -- The upper row index in the range of rows to print
          """
          for n in range(row_min, row_max + 1):
               for k in range(0, n + 1): 
                    sys.stdout.write("% 10d" % self.get_data(n, k))
               ## 
               print ""
          ##
     ## 

     def print_table(self):
          r""" 
          Prints the entire table of trangle values stored 
          locally by the instance of the object. 
          """
          self.print_table_rows(0, self.stored_rows)
     ##

## TriangularSequence 

r"""
Defines the triangle of unsigned Stirling numbers of the first kind (S1): 
"""
S1Triangle = TriangularSequence(1, 0, -1, 0, 0, 1)
def S1(n, k): return S1Triangle.get_data(n, k)

r"""
Defines the triangle of Stirling numbers of the second kind (S2): 
"""
S2Triangle = TriangularSequence(0, 1, 0, 0, 0, 1)
def S2(n, k): return S2Triangle.get_data(n, k)

r""" 
Defines the first-order Eulerian number triangle (E1): 
"""
E1Triangle = TriangularSequence(0, 1, 1, 1, -1, 0)
def E1(n, k): return E1Triangle.get_data(n, k)

r""" 
Defines the second-order Eulerian number triangle (E2): 
"""
E2Triangle = TriangularSequence(0, 1, 1, 2, -1, -1)
def E2(n, k): return E2Triangle.get_data(n, k)

r"""
Recursively defines the binomial coefficients (Binom), and 
a few useful variants of the symmetric index binomial coefficients, 
$\binom{n+k}{k}$, (BinomSymmetric) and the squared binomial 
coefficients (Binom2): 
"""
BinomTriangle = TriangularSequence(0, 0, 1, 0, 0, 1)
def Binom(n, k): return BinomTriangle.get_data(n, k)
def BinomSymmetric(n, k): return BinomTriangle.get_data(n + k, k)
def Binom2(n, k): return BinomTriangle.get_data(n, k) ** 2

def generating_function_coeff(ogf, z, n): 
     r""" 
     Returns the coefficient of $z^n$ in the Taylor series expansion 
     of the input ordinary generating function. 

     INPUTS: 
     - ``ogf`` -- The input ordinary generating function 
     - ``z`` -- The series variable in the power series expansions of ogf
     - ``n`` -- The nth coefficient index for some n >= 0
     """
     return taylor(ogf, z, 0, n + 1).coefficient(z, n)
##

def BernoulliB(n): 
     r""" 
     Implements another function generating the Bernoulli numbers, 
     $B_n$, by taking the coefficients of the sequence's exponential 
     generating function. 

     INPUTS: 
     - ``n`` -- The sequence index n >= 0 
     """
     z = var('z')
     bernoulli_egf = z / (exp(z) - 1)
     return factorial(n) * generating_function_coeff(bernoulli_egf, z, n)
## 

r"""
Define special SequenceGenerator objects for the guessing routines: 
"""
SeqGen_Bernoulli = \
     SequenceGenerator("NormalizedBernoulli", 
                       lambda n: factorial(n + 1) * BernoulliB(n), 
                       domain_dim = 1, latex_fmt = "({0}+1)! B_{{0}}")
SeqGen_FirstOrderHarmonic = \
     SequenceGenerator("NormalizedH1Number", 
                       lambda n: factorial(n) * harmonic_number(n, 1), 
                       domain_dim = 1, 
                       latex_fmt = "\\left({0}\\right)! H_{{0}}")
SeqGen_StirlingS1 = \
     SequenceGenerator("StirlingS1", 
                       lambda n, k: S1(n, k), 
                       domain_dim = 2, 
                       latex_fmt = "|s\\left({0}, {1}\\right)|")
SeqGen_StirlingS2 = \
     SequenceGenerator("StirlingS2", 
                       lambda n, k: S2(n, k), 
                       domain_dim = 2, 
                       latex_fmt = "S\\left({0}, {1}\\right)")
SeqGen_EulerianE1 = \
     SequenceGenerator("EulerianE1", 
                       lambda n, k: E1(n, k), 
                       domain_dim = 2, 
                       latex_fmt = "E_1\\left({0}, {1}\\right)")
SeqGen_EulerianE2 = \
     SequenceGenerator("EulerianE2", 
                       lambda n, k: E2(n, k), 
                       domain_dim = 2, 
                       latex_fmt = "E_2\\left({0}, {1}\\right)")
SeqGen_Binomial = \
     SequenceGenerator("Binom", 
                       lambda n, k: Binom(n, k), 
                       domain_dim = 2, 
                       latex_fmt = "\\binom{{0}}{{1}}")
SeqGen_BinomialSym = \
     SequenceGenerator("BinomSymmetric", 
                       lambda n, k: BinomSymmetric(n, k), 
                       domain_dim = 2, 
                       latex_fmt = "\\binom{{0}+{1}}{{1}}")
SeqGen_BinomialSquared = \
     SequenceGenerator("Binom2", 
                       lambda n, k: Binom2(n, k), 
                       domain_dim = 2, 
                       latex_fmt = "\\binom{{0}}{{1}}^2")
SeqGen_Factorial2 = \
     SequenceGenerator("Factorial2", 
                       lambda n: (2 ** n) * pochhammer(0.5, n) if (n % 2) == 1\
                                 else (2 ** n) * factorial(n), 
                       domain_dim = 1, 
                       latex_fmt = "\\left({0}\\right)!!")
SeqGen_NormalizedHarmonicNumber = \
     SequenceGenerator("NormalizedHNum", 
                       lambda n, r: factorial(n) * harmonic_number(n, r), 
                       domain_dim = 2, 
                       latex_fmt = "({0})! \\cdot H_{{0}}^{({1})}")

r"""
Default listing of the special SequenceGenerator objects available to the 
sequence formula guessing routines in Guess.py: 
""" 
SPECIAL_SEQUENCES = [ 
     SeqGen_Bernoulli, 
     SeqGen_FirstOrderHarmonic, 
     SeqGen_StirlingS1, 
     SeqGen_StirlingS2, 
     SeqGen_EulerianE1, 
     SeqGen_EulerianE2, 
     SeqGen_Binomial, 
     SeqGen_BinomialSym, 
     SeqGen_BinomialSquared, 
     SeqGen_Factorial2, 
     SeqGen_NormalizedHarmonicNumber, 
]

