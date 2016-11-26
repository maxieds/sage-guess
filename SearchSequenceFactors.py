#### SearchSequenceFactors.pyx : 
#### Implements a data structure for factorization-based searches of 
#### special factors of an (a sequence of) integer (integers)
#### Author:  Maxie D. Schmidt
#### Created: 2016.11.23 

from sage.all import *
from SequenceElement import *
from SpecialSequences import *
from GuessConfig import *

class SearchSequenceFactors(object):

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
               self.prime2index_functor = PrimeToIndexFunctor(num_primes)
          ##
     ## 

     def __len__(self): 
          return len(self._spseq_factors_list)
     ## 

     @property
     def num_primes(self): 
          return self._num_primes
     ##

     @num_primes.setter
     def num_primes(self, num_primes): 
          self._num_primes = num_primes
     ##

     @property
     def num_rows(self):
          return self._num_rows
     ## 

     @num_rows.setter
     def num_rows(self, num_rows):
          self._num_rows = num_rows
     ##

     @property
     def sequence_generator(self): 
          return self._sequence_generator
     ##

     @sequence_generator.setter
     def sequence_generator(self, seq_gen): 
          self._sequence_generator = seq_gen
     ## 

     def build_sequence_factors_list(self): 
          seqgen = self._sequence_generator
          seqgen.iter_first()
          seidx, sevalue = seqgen.next(self.num_rows)
          while sevalue != None:
               seq_elem = SequenceElement(sevalue, index_tuple = seidx)
               self._spseq_factors_list += [seq_elem]
               seidx, sevalue = seqgen.next(self.num_rows)
          ##
     ##

     def _compute_special_sequence_factors(self, seq_value, 
                                           map_index_pairs = False): 
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



