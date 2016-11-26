#### SequenceElement.py : 
#### Implements sequence elements, sequence generators, and related 
#### utility classes
#### Author:  Maxie D. Schmidt
#### Created: 2016.11.24

from sage.all import *
from GuessConfig import *

def AssertVariableType(var, vtype, exception_type = ValueError): 
     if not isinstance(var, vtype): 
          raise exception_type('Variable type must be ' + str(vtype) + '.')
     ##
##

class PrimePower(object): 
     
     def __init__(self, prime, ppow): 
          self._prime = prime
          self._prime_power = ppow
     ## 

     @property
     def prime(self): 
          return self._prime
     ##

     @property
     def ppow(self): 
          return self._prime_power
     ## 

## PrimePower

class PrimeToIndexFunctor(object): 

     def __init__(self, num_primes): 
          self.num_primes = num_primes
          self.prime2index_dict = \
               PrimeToIndexFunctor.build_prime2index_dict(num_primes)
     ## 
     
     def __len__(self): 
          return self.num_primes
     ## 

     def __call__(self, prime): 
          return self.prime2index_func(prime)
     ## 

     def __getitem__(self, prime): 
          return self.prime2index_func(prime)
     ## 

     def prime2index_func(self, prime): 
          if prime not in self.prime2index_dict:
               return -1
          else:
               return self.prime2index_dict[prime]
     ##

     @staticmethod
     def build_prime2index_dict(num_primes = DEFAULT_NUM_PRIMES): 
          prime2index_dict = dict()
          primes_list = primes_first_n(num_primes)
          for (pidx, prime) in enumerate(primes_list):
               prime2index_dict[prime] = pidx
          ## 
          return prime2index_dict
     ## 

     @staticmethod
     def index2prime_func(idx): 
          primes_list = primes_first_n(idx + 1)
          return primes_list[-1]
     ##

## PrimeToIndexFunctor

class LinearCombination(object): 

     def __init__(self, a, b): 
          self._a = a
          self._b = b
     ##

     @property
     def a(self): 
          return self._a
     ##

     @property
     def b(self): 
          return self._b
     ##

     def get_functor(self): 
          return lambda x: self.a * x + self.b
     ## 

     @staticmethod
     def get_functor(a, b): 
          return lambda x: a * x + b
     ## 

     def __call__(self, x): 
          return self.a * x + self.b
     ## 

     def __eq__(self, rhs): 
          if rhs == None: return False
          AssertVariableType(rhs, LinearCombination)
          return self.a == rhs.a and self.b == rhs.b
     ## 

     def __ne__(self, rhs): 
          if rhs == None: return True
          AssertVariableType(rhs, LinearCombination)
          return self.a != rhs.a or self.b != rhs.b
     ## 
     
     def __str__(self): 
          plus_minus_str = " + " if self.b > 0 else " - "
          bstr = "" if self.b == 0 else plus_minus_str + str(self.b)
          astr = "%s * " % str(self.a) if not IsOne(self.a) and \
                                          not IsZero(self.a) else ""
          rstr = "%s%s%s" % (astr, "%s", bstr)
          return rstr
     ## 

     def __latex__(self): 
          indetx_str = "\\left(%s\\right)"
          return str(self) % indetx_str
     ## 

     @staticmethod
     def find_linear_combination(index_seq, XY = X, check_all_data = True): 
          init_xpos = 0
          if init_xpos == len(index_seq) or init_xpos + 1 == len(index_seq):
               raise ValueError('Insufficient number of sequence elements.')
               return None
          ##
          datax, dataxp1 = index_seq[init_xpos][XY], index_seq[init_xpos + 1][XY]
          local_a = dataxp1 - datax
          local_b = datax - local_a * init_xpos
          lc = LinearCombination(local_a, local_b)
          if check_all_data:
               for (x, s) in enumerate(index_seq): 
                    if lc(x) != s: return None
               ##
          ## 
          return lc
     ##

     @staticmethod
     def find_all_linear_combinations(index_seq, domain_dim = 1, 
                                      check_all_data = True): 
          (E0Seq, E1Seq) = map(FiniteEnumeratedSet, index_seq[0:2])
          cprod = cartesian_product([E0Seq, E1Seq])
          lc_list = []
          for (idx0, idx1) in cprod: 
               lcx = LinearCombination.find_linear_combination([idx0, idx1], X, False)
               lcy = LinearCombination.find_linear_combination([idx0, idx1], Y, False)
               if lcx == None or lcy == None: continue
               lc_list += [(lcx, lcy)]
          ## 
          checked_lclist = []
          for lc in lc_list: 
               for sidx in range(2, len(index_seq)): 
                    lc_xindex, lc_yindex = lc[X](sidx), lc[Y](sidx)
                    if lc_xindex not in map(lambda t: t[X], index_seq[sidx]) or \
                       domain_dim == 2 and \
                       lc_yindex not in map(lambda t: t[Y], index_seq[sidx]): 
                         break
                    ##
               ## 
               checked_lclist += [lc]
          ## 
          return checked_lclist
     ## 

## LinearCombination

class SequenceElement(object): 
     
     def __init__(self, integer_value, index_tuple = None, 
                  num_primes = DEFAULT_NUM_PRIMES): 
          integer_value = Integer(integer_value)
          self._value = integer_value
          self._index_tuple = index_tuple
          self.num_primes = num_primes
          self.factor_result = [] if self.value == 0 else factor(integer_value)
          self._factor_list = list(self.factor_result)
          self.prime2index_func = PrimeToIndexFunctor(num_primes)
          self.ppow_hash_table = self.build_prime_power_hash(num_primes)
     ## 

     @property
     def value(self):
          return self._value
     ## 

     @property
     def index_tuple(self): 
          return self._index_tuple
     ##

     def __str__(self): 
          return str(self.value)
     ## 

     def __is_not__(self, rhs): 
          AssertVariableType(rhs, SequenceElement)
          return self.value != rhs.value
     ## 

     def __is__(self, rhs): 
          AssertVariableType(rhs, SequenceElement)
          return self.value != rhs.value
     ## 

     def __lt__(self, rhs): 
          AssertVariableType(rhs, SequenceElement)
          return self.value < rhs.value
     ## 

     def __le__(self, rhs): 
          AssertVariableType(rhs, SequenceElement)
          return self.value <= rhs.value
     ## 

     def __eq__(self, rhs): 
          AssertVariableType(rhs, SequenceElement)
          return self.value == rhs.value
     ## 

     def __ne__(self, rhs): 
          AssertVariableType(rhs, SequenceElement)
          return self.value != rhs.value
     ## 

     def __ge__(self, rhs): 
          AssertVariableType(rhs, SequenceElement)
          return self.value >= rhs.value
     ## 

     def __gt__(self, rhs): 
          AssertVariableType(rhs, SequenceElement)
          return self.value > rhs.value
     ## 

     def __contains__(self, ppower): 
          AssertVariableType(ppower, PrimePower)
          prime_index = self.prime2index_func(ppower.prime)
          if self.ppow_hash_table[prime_index] >= ppower.ppow: 
               return True
          else: 
               return False
          ##
     ##

     @staticmethod
     def div(a, b): 
          AssertVariableType(a, SequenceElement)
          AssertVariableType(b, SequenceElement)
          if IsZero(b.value) and IsZero(a.value): 
               return True
          elif IsZero(b.value): 
               return False
          elif IsOne(b.value) and IsOne(a.value):
               return True
          elif IsOne(b.value): 
               return False
          ## 
          for (prime, ppow) in b.factor_list():
               if PrimePower(prime, ppow) not in a: 
                    return False
               ##
          ##
          return True
     ## 

     def __idiv__(self, rhsb): 
          if not SequenceElement.div(self, rhsb): 
               raise RuntimeError('Division operation not an integer.')
               return None
          ## 
          div_value = self.value / rhsb.value
          return SequenceElement(div_value, None, self.num_primes)
     ## 

     def factor_list(self): 
          return self._factor_list
     ## 

     def unit(self): 
          return 0 if self.value == 0 else self.factor_result.unit()
     ## 

     def is_one(self): 
          return IsOne(abs(self.value))
     ## 

     def is_zero(self): 
          return IsZero(self.value)
     ## 

     def transform(self, tf_func): 
          tf_value = tf_func(self.value)
          self._value = tf_value
          self.ppow_hash_table = build_prime_power_hash(self.num_primes)
     ## 

     def build_prime_power_hash(self, num_primes): 
          ppow_hash = [0 for pidx in range(0, num_primes)]
          for (prime, ppow) in self.factor_list(): 
               prime_index = self.prime2index_func(prime)
               ppow_hash[prime_index] = ppow
          ## 
          return ppow_hash
     ## 

## SequenceElement

class SequenceGenerator(object): 
     
     def __init__(self, name, generator_func, domain_dim, latex_fmt = None): 
          self.name = name
          self.latex_fmt = latex_fmt
          self.functor = generator_func
          self.domain_dim = domain_dim
          self.current_index_tuple = (0, 0)
     ## 

     def get_print_string(self, fname): 
          var_inputs = ""
          if self.domain_dim == 1: 
               var_inputs = ('%s')
          elif self.domain_dim == 2:
               var_inputs = ('%s', '%s')
          else:
               raise IndexError('Unsupported domain dimension (should be 1 or 2).')
               return ""
          ##
          return fname % var_inputs 
     ## 

     def __str__(self): 
          if self.domain_dim == 1: 
               return self.get_print_string(self.name + '(%s)')
          elif self.domain_dim == 2: 
               return self.get_print_string(self.name + '(%s, %s)')
          else: 
               raise IndexError
               return ""
          ##
     ##

     def __latex__(self): 
          return self.get_print_string(self.latex_fmt)
     ## 

     def to_functor(self): 
          return self.functor
     ## 

     def __call__(self, index_tuple): 
          if index_tuple[0] < 0 or \
             self.domain_dim == 2 and index_tuple[1] < 0: 
               return Integer(0)
          ## 
          if self.domain_dim == 1: 
               return self.functor(index_tuple[0])
          elif self.domain_dim == 2: 
               return self.functor(index_tuple[0], index_tuple[1])
          else: 
               raise IndexError('Domain dimension should be 1 or 2.')
               return None
          ## 
     ## 
     
     def domain_dimension(self): 
          return self.domain_dim
     ##

     def Is1DFunc(self): return self.domain_dim == 1
     
     def Is2DFunc(self): return self.domain_dim == 2

     def iter_first(self): 
          self.current_index_tuple = (0, 0)
     ## 

     def next(self, num_rows): 
          if self.Is1DFunc(): 
               cur_x = self.current_index_tuple[0]
               next_x = self.current_index_tuple[0] + 1
               if next_x > num_rows: 
                    self.iter_first()
                    return self.current_index_tuple, None
               ## 
               self.current_index_tuple = (next_x, -1)
               return (cur_x, -1), self.functor(cur_x)
          elif self.Is2DFunc(): 
               n, k = self.current_index_tuple
               cur_n, cur_k = n, k
               if k < 0 or k >= n and n < num_rows: 
                    n += 1
               elif n < num_rows: 
                    k += 1
               else: 
                    self.iter_first()
                    return self.current_index_tuple, None
               self.current_index_tuple = (n, k)
               return (cur_n, cur_k), self.functor(cur_n, cur_k)
          else: 
               raise IndexError('Domain dimension should be 1 or 2.')
               return None
          ##
     ##

## SequenceGenerator 



