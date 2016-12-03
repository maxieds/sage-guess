
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
     product_factor_func = lambda i: x + i
     prod_func = prod(map(product_factor_func, range(0, n)))
     return prod(factor_list)
## def

def harmonic_number(n, r = 1):
     k = var('k')
     return sum(k ** -r, k, 1, n)
## def

def identity_func(): return lambda n: n
def zero_func(): return lambda n: 0
def one_func(): return lambda n, k = 1: 1
def constant_func(constant): return lambda n: constant
def linear_func(a, b): return lambda n: a * n + b
def sequence_func(seq_data):
     return lambda n: seq_data[n] if n >= 0 and n < len(seq_data) else 0

class triangular_sequence(object):

     def __init__(self, a, b, g, ap, bp, gp, max_rows = DEFAULT_NUM_ROWS):
          self.alpha = a 
          self.beta = b
          self.gamma = g
          self.alpha_prime = ap
          self.beta_prime = bp
          self.gamma_prime = gp
          self.stored_rows = max_rows
          self.rec_data = [[1]]
          self.generate_rows(1, max_rows)
     ## def 

     def generate_rows(self, row_min, row_max):
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
     ## def 

     def get_data(self, n, k): 
          if n < 0 or k < 0 or k > n or n > self.stored_rows:
               return 0
          elif n == 0 and k == 0: 
               return 1
          elif n > self.stored_rows:
               self.generate_rows(self.stored_rows + 1, n)
          ##
          return self.rec_data[n][k]
     ## def 

     def print_table_rows(self, row_min, row_max):
          #max_element = np.max(self.rec_data)
          for n in range(row_min, row_max + 1):
               for k in range(0, n + 1): 
                    sys.stdout.write("% 10d" % self.get_data(n, k))
               ## 
               print ""
          ##
     ## def 

     def print_table(self):
          self.print_table_rows(0, self.stored_rows)
     ## def

## class 
     
S1Triangle = triangular_sequence(1, 0, -1, 0, 0, 1)
S2Triangle = triangular_sequence(0, 1, 0, 0, 0, 1)
E1Triangle = triangular_sequence(0, 1, 1, 1, -1, 0)
E2Triangle = triangular_sequence(0, 1, 1, 2, -1, -1)
BinomTriangle = triangular_sequence(0, 0, 1, 0, 0, 1)

def S1(n, k): return S1Triangle.get_data(n, k)
def S2(n, k): return S2Triangle.get_data(n, k)
def E1(n, k): return E1Triangle.get_data(n, k)
def E2(n, k): return E2Triangle.get_data(n, k)
def Binom(n, k): return BinomTriangle.get_data(n, k)
def BinomSymmetric(n, k): return BinomTriangle.get_data(n + k, k)
def Binom2(n, k): return BinomTriangle.get_data(n, k) ** 2

def generating_function_coeff(ogf, z, n): 
     return taylor(ogf, z, 0, n + 1).coefficient(z, n)
##

def BernoulliB(n): 
     z = var('z')
     bernoulli_egf = z / (exp(z) - 1)
     return factorial(n) * generating_function_coeff(bernoulli_egf, z, n)
## 

#### Define special SequenceGenerator objects for the guessing routines: 
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

