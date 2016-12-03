
r"""
SequenceElement.py 

Implements sequence elements, sequence generators, and related utility classes. 

AUTHORS: 
- Maxie D. Schmidt (Created: 2016.11.24) 

"""

from sage.all import *
from GuessConfig import *

def AssertVariableType(var, vtype, exception_type = ValueError): 
     r""" 
     Assert function to mandate a variable type within the source code. 

     INPUTS: 
     - ``var`` -- A variable instance whose type is to be checked. 
     - ``vtype`` -- The required variable type. 
     - ``exception_type`` -- The exception thrown if the type of 
                             var does not match vtype. 
     """
     if not isinstance(var, vtype): 
          raise exception_type('Variable type must be ' + str(vtype) + '.')
     ##
##

class PrimePower(object): 
     r"""
     Class that represents a power of a prime. 
     """ 

     def __init__(self, prime, ppow): 
          r"""
          Initializes the PrimePower class instance. 

          INPUTS: 
          - ``prime`` -- An integer-typed prime number. 
          - ``ppow`` -- A non-negative integer power of the prime. 
          """ 
          self._prime = prime
          self._prime_power = ppow
     ## 

     def __int__(self): 
          r"""
          Defines the int(PrimePower object) type casting operator. 
          """
          return self.prime ** self.ppow
     ## 

     def __Integer__(self): 
          r"""
          Defines the Integer(PrimePower object) type casting operator. 
          """
          return Integer(self.prime ** self.ppow)
     ## 

     @property
     def prime(self): 
          r""" 
          Property representing the prime in the prime power object. 
          """ 
          return self._prime
     ##

     @property
     def ppow(self): 
          r"""
          Property representing the power in the prime power object. 
          """
          return self._prime_power
     ## 

## PrimePower

class PrimeToIndexFunctor(object): 
     r"""
     PrimeToIndexFunctor : 
     A class with a call method that returns the index $n-1$ of the nth prime. 
     This is useful in the setup of our factorization-based search routines. 
     """

     def __init__(self, num_primes): 
          r"""
          Initializes the PrimeToIndexFunctor object and its stored internal data.

          INPUTS: 
          - ``num_primes`` -- The number of distinct primes to include in the 
                              factorization-based searches. Locally, this means 
                              that we are only going to store lookup data for the 
                              first num_primes prime numbers. 
          """
          self.num_primes = num_primes
          self.prime2index_dict = \
               PrimeToIndexFunctor.build_prime2index_dict(num_primes)
     ## 
     
     def __len__(self): 
          r""" 
          Returns the number of primes we are considering in the lookup methods. 
          """
          return self.num_primes
     ## 

     def __call__(self, prime): 
          r"""
          Returns the array index, $n-1$, for the nth prime denoted by prime, 
          or -1 if this prime is out of index in the local searches. 
          A convenient operator wrapper around the function 
          prime2index_func defined in the source below. 

          INPUTS: 
          - ``prime`` -- A prime number. 
          """
          return self.prime2index_func(prime)
     ## 

     def __getitem__(self, prime): 
          r"""
          Returns the array index, $n-1$, for the nth prime denoted by prime, 
          or -1 if this prime is out of index in the local searches. 
          A convenient operator wrapper around the function 
          prime2index_func defined in the source below. 

          INPUTS: 
          - ``prime`` -- A prime number. 
          """
          return self.prime2index_func(prime)
     ## 

     def prime2index_func(self, prime): 
          r"""
          Returns the array index, $n-1$, for the nth prime denoted by prime, 
          or -1 if this prime is out of index in the local searches. 

          INPUTS: 
          - ``prime`` -- A prime number. 
          """
          if prime not in self.prime2index_dict:
               return -1
          else:
               return self.prime2index_dict[prime]
     ##

     @staticmethod
     def build_prime2index_dict(num_primes = DEFAULT_NUM_PRIMES): 
          r""" 
          Builds the local prime to index mapping dictionary. 

          INPUTS: 
          - ``num_primes`` -- The number of distinct primes to include in the 
                              dictionary lookups. Inputs of primes larger than 
                              $p_{num_primes}$ cause the lookup functions / 
                              operators defined above to return -1. 
          """
          prime2index_dict = dict()
          primes_list = primes_first_n(num_primes)
          for (pidx, prime) in enumerate(primes_list):
               prime2index_dict[prime] = pidx
          ## 
          return prime2index_dict
     ## 

     @staticmethod
     def index2prime_func(idx): 
          r"""
          Given an index i, returns the ith prime number. 

          INPUTS: 
          - ``idx`` -- An integer index into the sequence of primes. 
          """ 
          primes_list = primes_first_n(idx + 1)
          return primes_list[-1]
     ##

## PrimeToIndexFunctor

class LinearCombination(object): 
     r"""
     LinearCombination : 
     Represents a linear index combination of the form ax+b for an input 
     index x where a, b are integer-valued scalars (not both zero). 
     """

     def __init__(self, a, b): 
          r"""
          Initializes the LinearCombination object. 

          INPUTS: 
          - ``a`` -- The leading coefficient of the combination ax+b
          - ``b`` -- The constant term of the combination ax+b
          """
          self._a = a
          self._b = b
     ##

     @property
     def a(self): 
          r"""
          Property returning the leading coefficient of the combination ax+b
          """
          return self._a
     ##

     @property
     def b(self): 
          r"""
          Property returning the constant coefficient of the combination ax+b
          """
          return self._b
     ##

     def get_functor(self): 
          r""" 
          Returns an unevaluated lambda function representing the 
          linear combination. 
          """
          return lambda x: self.a * x + self.b
     ## 

     @staticmethod
     def get_functor(a, b): 
          r""" 
          Returns an unevaluated lambda function representing the 
          linear combination ax+b for an indeterminate x. 
          
          INPUTS:
          - ``a`` -- The leading coefficient of the combination
          - ``b`` -- The constant term in the combination
          """
          return lambda x: a * x + b
     ## 

     def __call__(self, x): 
          r"""
          Evaluates the linear combination ax+b at x. 

          INPUTS: 
          - ``x` -- The input index to the combination function 
          """
          return self.a * x + self.b
     ## 

     def __eq__(self, rhs): 
          r""" 
          Tests whether the current linear combination is equal to the rhs input. 

          INPUTS: 
          - ``rhs`` -- None (representing a NULL combination object), or
                       another LinearCombination 
          """
          if rhs == None: return False
          AssertVariableType(rhs, LinearCombination)
          return self.a == rhs.a and self.b == rhs.b
     ## 

     def __ne__(self, rhs): 
          r""" 
          Tests whether the current linear combination is unequal to the rhs input. 

          INPUTS: 
          - ``rhs`` -- None (representing a NULL combination object), or
                       another LinearCombination 
          """
          if rhs == None: return True
          AssertVariableType(rhs, LinearCombination)
          return self.a != rhs.a or self.b != rhs.b
     ## 
     
     def __str__(self): 
          r""" 
          Returns a parameterized string representation of the object. 

          EXAMPLES: 
          
          The indeterminate x in ax+b is left unevaluated as a "%s%" 
          placeholder so that the string representation can in this 
          sense be "evaluated" at a generic, or explicit, value of x: 

          sage: from SequenceElement import LinearCombination
          sage: lc = LinearCombination(2, 1) # 2x+1
          sage: lcstr = str(lc)
          sage: nx = var('x')
          sage: print lcstr % nx
          2 * x + 1
          sage: print lcstr % 3 # should represent 7
          2 * 3 + 1
          
          Now we can create a function that returns the correctly 
          evaluated (or unevaluated as the case may be) string: 

          sage: from GuessConfig import IsInteger
          sage: eval_func = lambda nx: str(eval(lcstr % nx)) if IsInteger(nx) \
          ....:                        else lcstr % nx
          sage: eval_func(nx)
          '2 * x + 1'
          sage: eval_func(4) 
          '9' 

          """ 
          plus_minus_str = " + " if self.b > 0 else " - "
          bstr = "" if self.b == 0 else plus_minus_str + str(self.b)
          astr = "%s * " % str(self.a) if not IsOne(self.a) and \
                                          not IsZero(self.a) else ""
          rstr = "%s%s%s" % (astr, "%s", bstr)
          return rstr
     ## 

     @staticmethod
     def find_linear_combination(index_seq, XY = X, check_all_data = True): 
          r""" 
          Finds a linear combination of the index sequence of tuples indexed 
          either over the first component of the sequence (XY = X), or 
          over the second component of the sequence (XY = Y). 

          INPUTS: 
          - ``index_seq`` -- A list of 2-tuples 
          - ``XY`` -- Either X (find linear combination over the x components)
                      or Y (over the y components) 
          - ``check_all_data`` -- Whether to check that the returned linear 
                                  combination matches all of the input 
                                  index_seq elements, or just the first two. 
                                  Only checking the first two values 
                                  (option = False) allows all of the possible 
                                  linear combinations generated by 
                                  find_all_linear_combinations below to be 
                                  computed quickly and left to verify later. 
          """ 
          init_xpos = 0
          if init_xpos == len(index_seq) or init_xpos + 1 == len(index_seq):
               raise ValueError('Insufficient number of sequence elements.')
               return None
          ##
          datax, dataxp1 = index_seq[init_xpos][XY], index_seq[init_xpos + 1][XY]
          local_a = dataxp1 - datax
          if local_a < 0: return None
          local_b = datax - local_a * init_xpos
          lc = LinearCombination(local_a, local_b)
          if check_all_data:
               for (x, s) in enumerate(index_seq): 
                    if lc(x) != s[XY]: return None
               ##
          ## 
          return lc
     ##

     @staticmethod
     def find_all_linear_combinations(index_seq, domain_dim = 1): 
          r"""
          Finds all possible accurate index combinations formed by the 
          input sequence of index 2-tuples. 

          INPUTS: 
          - ``index_seq`` -- A list of lists of index 2-tuples 
          - ``domain_dim`` -- Whether the corresponding special sequence is 
                              indexed by the index_seq is one or two-dimensional. 
          """ 
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
          r""" 
          Initializes the SequenceElement object. 

          INPUTS: 
          - ``integer_value`` -- The value associated with this sequence element 
          - ``index_tuple`` -- The 2-tuple corresponding to the sequence index 
                               of this sequence element (if defined) 
          - ``num_primes `` -- The number of distinct primes to use in 
                               checking integer divisibility properties. 
          """ 
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
          r""" 
          Property returning the integer value of the sequence element. 
          """
          return self._value
     ## 

     @property
     def index_tuple(self): 
          r""" 
          Property returning the index-tuple associated with this 
          sequence element. 
          """
          return self._index_tuple
     ##

     def __str__(self): 
          r""" 
          Returns a string representation of the sequence element. 
          """
          return str(self.value)
     ## 

     def __int__(self): 
          r""" 
          Implements a type cast of the sequence element to a Python int. 
          """
          return int(self.value)
     ##

     def __Integer__(self): 
          r""" 
          Implements a type cast of the sequence element to a Sage Integer. 
          """
          return Integer(self.value)
     ##

     def __is_not__(self, rhs): 
          r""" 
          Tests whether the local object is not equal to the RHS input. 

          INPUTS: 
          - ``rhs`` -- Right-hand-side operand of type SequenceElement 
          """
          AssertVariableType(rhs, SequenceElement)
          return self.value != rhs.value
     ## 

     def __is__(self, rhs): 
          r""" 
          Tests whether the local object is equal to the RHS input. 

          INPUTS: 
          - ``rhs`` -- Right-hand-side operand of type SequenceElement 
          """
          AssertVariableType(rhs, SequenceElement)
          return self.value != rhs.value
     ## 

     def __lt__(self, rhs): 
          r""" 
          Tests whether the local object value is less than to the RHS input. 

          INPUTS: 
          - ``rhs`` -- Right-hand-side operand of type SequenceElement 
          """
          AssertVariableType(rhs, SequenceElement)
          return self.value < rhs.value
     ## 

     def __le__(self, rhs): 
          r""" 
          Tests whether the local object value is less than or 
          equal to the RHS input. 

          INPUTS: 
          - ``rhs`` -- Right-hand-side operand of type SequenceElement 
          """
          AssertVariableType(rhs, SequenceElement)
          return self.value <= rhs.value
     ## 

     def __eq__(self, rhs): 
          r""" 
          Tests whether the local object value is equal to the RHS input. 

          INPUTS: 
          - ``rhs`` -- Right-hand-side operand of type SequenceElement 
          """
          AssertVariableType(rhs, SequenceElement)
          return self.value == rhs.value
     ## 

     def __ne__(self, rhs): 
          r""" 
          Tests whether the local object value is not equal to the RHS input. 

          INPUTS: 
          - ``rhs`` -- Right-hand-side operand of type SequenceElement 
          """
          AssertVariableType(rhs, SequenceElement)
          return self.value != rhs.value
     ## 

     def __ge__(self, rhs): 
          r""" 
          Tests whether the local object value is greater than or 
          equal to the RHS input. 

          INPUTS: 
          - ``rhs`` -- Right-hand-side operand of type SequenceElement 
          """
          AssertVariableType(rhs, SequenceElement)
          return self.value >= rhs.value
     ## 

     def __gt__(self, rhs): 
          r""" 
          Tests whether the local object value is greater than to the RHS input. 

          INPUTS: 
          - ``rhs`` -- Right-hand-side operand of type SequenceElement 
          """
          AssertVariableType(rhs, SequenceElement)
          return self.value > rhs.value
     ## 

     def __contains__(self, ppower): 
          r""" 
          Operator returns a boolean indicator of whether the input 
          prime power divides the local object value. 

          INPUTS: 
          - ``ppower`` -- An input initialized PrimePower object 
          """ 
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
          r""" 
          Tests whether the integer a is divisible by the integer b 
          with some sensible special cases implementing the logic of 
          returning intended special sequence factors. 

          INPUTS: 
          - ``a`` -- An integer-valued SequenceElement 
          - ``b`` -- An integer-valued SequenceElement 
          """ 
          AssertVariableType(a, SequenceElement)
          AssertVariableType(b, SequenceElement)
          if IsZero(b.value) and IsZero(a.value): 
               return True
          elif IsZero(b.value): 
               return False
          elif IsOne(b.value) and IsOne(a.value):
               return True
          elif IsOne(b.value): 
               return True
          ## 
          for (prime, ppow) in b.factor_list():
               if PrimePower(prime, ppow) not in a: 
                    return False
               ##
          ##
          return True
     ## 

     def __idiv__(self, rhsb): 
          r""" 
          Implements the /= operator where division is tested by the 
          function div(a, b) in the source above. 

          INPUTS: 
          - ``rhsb`` -- A right-hand-side SequenceElement operand
          """
          if not SequenceElement.div(self, rhsb): 
               raise RuntimeError('Division operation not an integer.')
               return None
          ## 
          div_value = self.value / rhsb.value
          return SequenceElement(div_value, None, self.num_primes)
     ## 

     def factor_list(self): 
          r""" 
          Returns a list of factors of the value associated with this 
          object in the form of (prime, prime_power). 
          """
          return self._factor_list
     ## 

     def unit(self): 
          r"""
          Tests whether the value associated with this object is a unit. 
          """
          return 0 if self.value == 0 else self.factor_result.unit()
     ## 

     def is_one(self): 
          r""" 
          Tests whether the absolute value of the value associated with 
          this object is equal to one. 
          """
          return IsOne(abs(self.value))
     ## 

     def is_zero(self): 
          r"""
          Tests whether the value associated with this object is zero. 
          """
          return IsZero(self.value)
     ## 

     def transform(self, tf_func): 
          r""" 
          Transforms the value of the object according to the input 
          transformation function. 
          Eventually, we will implement more sequence transformations, 
          including the generating function transformation operations 
          implemented in the GFun.py Sage package. 

          INPUTS: 
          - ``tf_func`` -- A transformation function returning an integer value 
          """
          tf_value = tf_func(self.value)
          self._value = tf_value
          self.ppow_hash_table = build_prime_power_hash(self.num_primes)
     ## 

     def build_prime_power_hash(self, num_primes): 
          r""" 
          Internal helper function that 
          constructs an array mapping primes (indexed by their sequence order) 
          to non-negative integer powers of these primes. 

          INPUTS: 
          - ``num_primes`` -- The number of distinct primes to consider in the 
                              hash mapping constructed here. 
                              This value should be set so that the largest 
                              possible prime sequence divisor is less than or 
                              equal to $p_{num_primes}$. 
          """
          ppow_hash = [0 for pidx in range(0, num_primes)]
          for (prime, ppow) in self.factor_list(): 
               prime_index = self.prime2index_func(prime)
               ppow_hash[prime_index] = ppow
          ## 
          return ppow_hash
     ## 

## SequenceElement

class SequenceGenerator(object): 
     r""" 
     SequenceGenerator : 
     Class that defines and implements calls to the values of a special 
     sequence. The stock SPECIAL_SEQUENCES list defined in 
     SpecialSequences.py provides examples of the usage of these objects in 
     practice with our sequence guessing routines. 
     """

     def __init__(self, name, generator_func, domain_dim, latex_fmt = None): 
          r""" 
          Initializes the SequenceGenerator object. 

          INPUTS: 
          - ``name`` -- The print name identifier of the sequence 
          - ``generator_func`` -- A one or two-argument function that 
                                  generates the special sequence over its
                                  domain 
          - ``domain_dim`` -- The dimension of the domain of the sequence. 
                              Only 1 and 2 are supported. 
          - ``latex_fmt`` -- A format string giving a LaTeX representation 
                             of the sequence 
          """
          self.name = name
          self.latex_fmt = latex_fmt
          self.functor = generator_func
          self.domain_dim = domain_dim
          self.current_index_tuple = (0, 0)
     ## 

     def get_print_string(self, fname): 
          r"""
          Returns a string representation of the special sequence 
          represented by this object. The inputs are left "unevaluated" as 
          "%s" placeholders for the actual inputs to the sequence to be 
          printed by external functions. 
          """
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
          r"""
          Returns a string representation of the special sequence 
          represented by this object. The inputs are left "unevaluated" as 
          "%s" placeholders for the actual inputs to the sequence to be 
          printed by external functions. 
          """
          if self.domain_dim == 1: 
               return self.get_print_string(self.name + '(%s)')
          elif self.domain_dim == 2: 
               return self.get_print_string(self.name + '(%s, %s)')
          else: 
               raise IndexError
               return ""
          ##
     ##

     def to_functor(self): 
          r""" 
          Returns the function generating the sequence. 
          Typically, but not always, defined as an unevaluated lambda function. 
          """
          return self.functor
     ## 

     def __call__(self, index_tuple): 
          r"""
          Returns the value of the special sequence at the specified index inputs. 

          INPUTS: 
          - ``index_tuple`` -- A 2-tuple representing the one-dimensional 
                               (in which case the second tuple input is ignored), 
                               or the two-dimensional sequence index inputs 
          """
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
          r"""
          Returns the dimension of the domain of the special sequence 
          represented by this object. Currently only one and two-dimensional 
          sequences are supported by this class, and by the sequence formula 
          guessing routines in the package. 
          """
          return self.domain_dim
     ##

     def Is1DFunc(self): 
          r""" 
          Tests whether the special sequence represented by the object is 
          one-dimensional
          """
          return self.domain_dim == 1
     ##
     
     def Is2DFunc(self): 
          r"""
          Tests whether the special sequence represented by the object is 
          two-dimensional. 
          """
          return self.domain_dim == 2
     ## 

     def iter_first(self): 
          r""" 
          Sets the internal iterator over the sequence values to its 
          initial position. 
          """
          self.current_index_tuple = (0, 0)
     ## 

     def next(self, num_rows): 
          r""" 
          Advances the internal sequence value iterator and returns a pair of the 
          current sequence index 2-tuple and the sequence value at this position. 
          When the maximum number of sequence elements specified by the input 
          is exceeded, the value in the returned pair is given by None. 

          INPUTS: 
          - ``num_rows`` -- The number of rows for two-dimensional triangular 
                            sequences, or alternately, the maximum number of 
                            sequence elements in a one-dimensional sequence. 

          EXAMPLES: 

          Obtain a list of the elements of the first three rows of the 
          unsigned triangle of Stirling numbers of the first kind: 

          sage: from SpecialSequences import SeqGen_StirlingS1
          sage: max_rows = 4; s1seqgen = SeqGen_StirlingS1
          sage: s1seqgen.iter_first()
          sage: seq_elems = []
          sage: while True: 
          ....:      idx_tuple, seq_value = s1seqgen.next(max_rows)
          ....:      if seq_value == None: break
          ....:      seq_elems += [seq_value] 
          ....:      print "S1%s = %d" % (idx_tuple, seq_value)
          S1(0, 0) = 1
          S1(1, 0) = 0
          S1(1, 1) = 1
          S1(2, 0) = 0
          S1(2, 1) = 1
          S1(2, 2) = 1
          S1(3, 0) = 0
          S1(3, 1) = 2
          S1(3, 2) = 3
          S1(3, 3) = 1
          sage: print seq_elems
          [1, 0, 1, 0, 1, 1, 0, 2, 3, 1]
          
          """
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
                    k = 0
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



