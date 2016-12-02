
r"""
SearchSequenceFactors.py 

Implements a data structure for factorization-based searches of 
special factors of an (a sequence of) integer (integers)

AUTHORS: 
- Maxie D. Schmidt (Created: 2016.11.23)

"""

from sage.all import *
from SequenceElement import *
from SpecialSequences import *
from GuessConfig import *

class SearchSequenceFactors(object):
     r"""
     SearchSequenceFactors 

     Class that implements searching of an input sequence for 
     special sequence factors corresponding to a SequenceGenerator 
     defined in the source of SpecialSequences.py. 
     """

     def __init__(self, num_primes = DEFAULT_NUM_PRIMES, 
                  num_rows = DEFAULT_NUM_ROWS, 
                  seqgen = None, build_immediately = True): 
          self._num_primes = num_primes
          self._num_rows = num_rows
          self._sequence_generator = seqgen
          self._spseq_factors_list = []
          self._prime2index_functor = None
          if build_immediately: 
               self.build_sequence_factors_list()
          ##
     ## 

     def __len__(self): 
          r""" 
          The number of special sequence factors in the search list. 
          """
          return len(self._spseq_factors_list)
     ## 

     @property
     def num_primes(self): 
          r"""
          Returns the number of primes used in the factorization-based 
          searches for special sequence factors. 
          """
          return self._num_primes
     ##

     @num_primes.setter
     def num_primes(self, num_primes): 
          r"""
          Sets the num_primes property of the class. 

          INPUTS: 
          - ``num_primes`` -- Set the number of primes to use in the 
                              factorization-based searches for special 
                              sequence factors. 
                              Class defaults to DEFAULT_NUM_PRIMES (=64). 
          """
          self._num_primes = num_primes
     ##

     @property
     def num_rows(self):
          r""" 
          Defines the number of sequence elements, or triangular 
          special sequence rows, used in the search for special 
          sequence factors. 
          """
          return self._num_rows
     ## 

     @num_rows.setter
     def num_rows(self, num_rows):
          r"""
          Sets num_rows property of the class. 

          INPUTS: 
          - ``num_rows`` -- Positive integer-valued setting for the 
                            number of special sequence rows (1D elements).
          """
          self._num_rows = num_rows
     ##

     @property
     def sequence_generator(self): 
          r"""
          Returns the SequenceGenerator object used to generate the 
          local special sequence factors. 
          """
          return self._sequence_generator
     ##

     @sequence_generator.setter
     def sequence_generator(self, seq_gen): 
          r"""
          Sets the sequence_generator property of the class. 

          INPUTS: 
          - ``seq_gen`` -- A special sequence represented by a 
                           SequenceGenerator object. 
          """
          self._sequence_generator = seq_gen
     ## 

     def build_sequence_factors_list(self): 
          r"""
          Builds the locally-stored list of special sequence factors by 
          iterating over the elements up to and including the last element 
          specified by self.num_rows. 
          
          EXAMPLES: 

          To manually setup the class instance, you may set the 
          num_primes, num_rows, and sequence_generator properties, and then 
          call this function: 
          
          sage: from SpecialSequences import SeqGen_StirlingS1
          sage: ssf = SearchSequenceFactors(build_immediately = False)
          sage: ssf.num_primes = 128
          sage: ssf.num_rows = 24
          sage: ssf.sequence_generator = SeqGen_StirlingS1
          sage: ssf.build_sequence_factors_list()
          
          """
          seqgen = self._sequence_generator
          seqgen.iter_first()
          seidx, sevalue = seqgen.next(self.num_rows)
          while sevalue != None:
               seq_elem = SequenceElement(sevalue, index_tuple = seidx)
               self._spseq_factors_list += [seq_elem]
               seidx, sevalue = seqgen.next(self.num_rows)
          ##
          self.prime2index_functor = PrimeToIndexFunctor(self.num_primes)
     ##

     def _compute_special_sequence_factors(self, seq_value, 
                                           map_index_pairs = False): 
          r""" 
          Computes a list of special sequence factors dividing the input 
          sequence value. 

          INPUTS: 
          - ``seq_value`` -- A particular (integer) value from the input sequence 
          - ``map_index_pairs`` -- Specifies whether to return a list of the 
                                   corresponding SequenceGenerator index 
                                   tuple for each special sequence factor 
                                   of seq_value found in the search. 
          """
          svalue = SequenceElement(seq_value)
          factors_list = self._spseq_factors_list
          rfactors_list = []
          for (ssidx, ssvalue) in enumerate(factors_list): 
               if SequenceElement.div(svalue, ssvalue): 
                    rfactors_list += [ssvalue]
               ##
          ##
          if map_index_pairs: 
               rfactors_list = map(lambda se: se.index_tuple, rfactors_list)
          ##
          return rfactors_list
     ## 

     def compute_remaining_terms(self, sequence, lcidx):
          r""" 
          Computes a list of remaining terms in the input sequence after 
          removing the special sequence factors indexed by a 
          LinearCombination 2-tuple. 

          INPUTS: 
          - ``sequence`` -- The original input sequence (list of integers) 
          - ``lcidx`` -- The special sequence indexing function returned by a 
                         call to LinearCombination.find_all_linear_combinations. 
          """
          rem_sequence = []
          for (x, selem) in enumerate(sequence): 
               sseq_value = self.sequence_generator((lcidx[X](x), lcidx[Y](x)))
               selem = 0 if IsZero(selem) or IsZero(sseq_value)\
                         else selem / sseq_value 
               rem_sequence += [selem]
          ##
          return rem_sequence
     ## 

     def compute_special_sequence_factors(self, sequence): 
          r""" 
          Computes a complete list of special sequence factors specified by the 
          local SequenceGenerator object and a 2-tuple of LinearCombination 
          index objects. Also, returns for each matching index tuple a list of the 
          remaining sequence terms after the special sequence factors are 
          removed from the input sequence elements. 

          INPUTS: 
          - ``sequence`` -- The original input sequence (list of integers) 
          """ 
          seq_factors_list = []
          for seq_elem in sequence: 
               cur_elem_factors = \
                    self._compute_special_sequence_factors(seq_elem, 
                                                           map_index_pairs = True)
               seq_factors_list += [cur_elem_factors]
          ## 
          lc_list = LinearCombination.find_all_linear_combinations(\
                    seq_factors_list, self.sequence_generator.domain_dimension())
          lc_list = map(lambda lc: \
                    (lc, self.compute_remaining_terms(sequence, lc)), lc_list)
          return lc_list
     ##

## SearchSequenceFactors 

