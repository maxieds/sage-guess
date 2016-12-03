
r""" 
Guess.py 

The main collaborative sequence formula guessing routines in Sage
(and some newly added routines to handle factors of special sequences)

AUTHORS: 
- Maxie D. Schmidt (Created: 2016.11.24)

"""

from sage.all import *
from GuessingFormulasForSequences import *
from SpecialSequences import *
from SequenceElement import AssertVariableType
from SearchSequenceFactors import *
from ColorizeString import print_color_string
from GuessConfig import *

class OEISFormula(GuessFormulaResult): 
     r"""
     OEISFormula : 
     Represents the match data from a matching sequence 
     found in the OEIS database using the function 
     guess_oeis_formula defined in the source below. 
     """

     def __init__(self, oeis_search_result): 
          r""" 
          Initializes the OEISFormula object. 

          INPUTS: 
          - ``oeis_search_result`` -- A sequence search result returned from
                                      the built-in Sage oeis(...) function
          """
          self._seq_number = oeis_search_result.id()
          self._name = oeis_search_result.name()
          self.functor = sequence_func(list(oeis_search_result.first_terms()))
     ## 

     @property
     def sequence_number(self): 
          r""" 
          Property returning the matching sequence's OEIS sequence 
          number, AXXXXXX. 
          """
          return self._seq_number
     ## 

     @property
     def name(self): 
          r""" 
          Property returning the name of the OEIS sequence match. 
          """
          return self._name
     ## 

     def __str__(self): 
          r""" 
          Returns a string representation of the OEIS sequence match. 
          """
          return "%s(n)" % self.sequence_number
     ## 

     def __eq__(self, rhs): 
          r""" 
          Tests whether the local OEIS sequence is equal to the 
          right-hand-side sequence by comparing sequence numbers. 
          
          INPUTS: 
          - ``rhs`` -- Right-hand-side OEISFormula operand 
          """
          AssertVariableType(rhs, OEISFormula)
          return self.sequence_number == rhs.sequence_number
     ## 

     def __ne__(self, rhs): 
          r""" 
          Tests whether the local OEIS sequence is not equal to the 
          right-hand-side sequence by comparing sequence numbers. 
          
          INPUTS: 
          - ``rhs`` -- Right-hand-side OEISFormula operand 
          """
          AssertVariableType(rhs, OEISFormula)
          return self.sequence_number != rhs.sequence_number
     ## 

     def desc(self): 
          r""" 
          Returns a description of the sequence match. 
          """
          return "%s: %s" % (self.sequence_number, self.name)
     ## 

     def print_summary(self): 
          r""" 
          Prints a summary of the OEIS sequence match data. 
          """
          print "   >>", self.desc()
          print "   >> First Sequence Values: ", \
                [self.functor(x) for x in range(0, 10)]
     ## 

## 

def guess_oeis_formula(sequence): 
     r""" 
     Returns a list of matches for the input sequence in the OEIS database. 

     INPUTS: 
     - ``sequence`` -- A list of integers 

     EXAMPLES: 

     Find matches in the OEIS database for the first few triangular numbers: 

     sage: seqnums = guess_oeis_formula([1, 3, 6, 10, 15, 21])
     sage: map(lambda sm: sm.print_summary(), seqnums)
        >> A000217: Triangular numbers: a(n) = binomial(n+1,2) = n(n+1)/2 = 0 + 1 + 2 + ... + n.
        >> First Sequence Values:  [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]
        >> A037123: a(n) = a(n-1) + Sum of digits of n.
        >> First Sequence Values:  [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]
     [None, None]

     """
     search_results = oeis(sequence, max_results = MAX_OEIS_RESULTS)
     search_results = map(OEISFormula, search_results)
     return search_results
## 

class GeneratingFunctionMatch(GuessFormulaResult): 
     r"""
     GeneratingFunctionMatch : 
     Represents a formula match in the form of a generating function 
     enumerating the sequence obtained with guess_generating_function 
     defined in the source below. 
     """

     def __init__(self, sequence, gfmatch): 
          r""" 
          Initializes the GeneratingFunctionMatch object. 

          INPUTS: 
          - ``sequence`` -- An input sequence as a list of values 
          - ``gfmatch`` -- A generating function enumerating the sequence 
                           obtained from a call to guess_generating_function 
          """
          self.gfmatch = gfmatch
          self.functor = sequence_func(sequence)
     ## 

     def __str__(self): 
          r"""
          Returns a string representation of the match object. 
          """
          return "[x^n] " + str(self.gfmatch.ogf())
     ## 
     
     def print_summary(self): 
          r""" 
          Prints a summary of the matching generating function. 
          """
          print "   >> Sequence OGF: ", self.gfmatch.ogf(), " = ", \
                self.gfmatch.series(4)
          print "   >> Recurrence: ", self.gfmatch.recurrence_repr()
     ## 


## GeneratingFunctionMatch

def guess_generating_function(sequence): 
     r""" 
     Guesses a formula for a generating function enumerating the 
     input sequence. 

     INPUTS: 
     - ``sequence`` -- A list of integer-like terms 
     
     EXAMPLES: 
     We can guess a generating function for particular columns of 
     special triangular sequences. In this example, we guess a generating 
     function for the third column of the Stirling numbers of the 
     second kind: 

     sage: from SpecialSequences import *
     sage: s2seq = [S2(n, 3) for n in range(3, 12)]
     sage: print s2seq
     [1, 6, 25, 90, 301, 966, 3025, 9330, 28501]
     sage: guess_generating_function(s2seq)[0].print_summary()
        >> Sequence OGF:  1/(-6*x^3 + 11*x^2 - 6*x + 1)  =  1 + 6*x + 25*x^2 + 90*x^3 + O(x^4)
        >> Recurrence:  Homogenous linear recurrence with constant coefficients of degree 3: a(n+3) = 6*a(n+2) - 11*a(n+1) + 6*a(n), starting a(0...) = [1, 6, 25]
     
     Compare with the corresponding result returned from a match in the 
     OEIS database: 

     sage: guess_oeis_formula(s2seq)[0].print_summary()
        >> A000392: Stirling numbers of second kind S(n,3).
        >> First Sequence Values:  [0, 0, 0, 1, 6, 25, 90, 301, 966, 3025]
     
     Note that the corresponding exponential generating function of the 
     sequence, despite having a well-known closed-form expression, is 
     not found so readily by these functions. 
     
     """
     C, x = CFiniteSequences(QQ, 'x').objgen()
     smatch = C.guess(sequence) 
     if smatch == 0: 
          return []
     else: 
          return [GeneratingFunctionMatch(sequence, smatch)]
     ## 
## 

def guess_remaining_sequence(sequence): 
     r""" 
     Guesses a sequence formula by defaulting to the subroutines we have 
     packaged as the GuessingFormulasForSeqs functions, the 
     guess_oeis_formula, and the guess_generating_function functions defined 
     above. 
     
     This function is designed to process remaining terms of a sequence 
     after removing the more difficult to guess special sequence factors 
     identified by the guess(...) function defined below. This function 
     should obtain closed-form formula matches for most already easily 
     identified "remaining" sequences. 

     INPUTS: 
     - ``sequence`` -- A list of integer-like values 
     """
     guess_pkg_results = GuessingFormulasForSeqs.guess(sequence)
     #guess_pkg_results += guess_oeis_formula(sequence)
     #guess_pkg_results += guess_generating_function(sequence)
     return guess_pkg_results 
## 

class GuessSpecialSeqResult(GuessFormulaResult): 
     r""" 
     A wrapper class packaging the data for a sequence formula match 
     containing a list of special function factors indexed over some 
     prescribed 2-tuple of linear combinations of the sequence index. 
     """ 

     def __init__(self, remseq, seqgen, lcidx_func): 
          r""" 
          Initializes the GuessSpecialSeqResult object with the data from 
          a first matching special sequence factor. 

          INPUTS: 
          - ``remseq`` -- The leftover remaining sequence after the identified 
                          special sequence factors have been removed from the 
                          input sequence 
          - ``seqgen`` -- A SequenceGenerator object generating the first 
                          special sequence factor 
          - ``lcidx_func`` -- A 2-tuple of LinearCombination objects 
                              representing the inputs to the first special 
                              sequence factor (in the list of all 
                              special sequence factors)
          """
          self.spseq_factors = [(seqgen, lcidx_func)]
          self.remseq = remseq
          self.finalize = False
     ## 

     @staticmethod
     def merge(match1, match2, combined_remseq): 
          r""" 
          Merges the special sequence factor data between two 
          GuessSpecialSeqResult objects and returns the combined result. 

          INPUTS: 
          - ``match1`` -- Initialized GuessSpecialSeqResult object 
          - ``match2`` -- Another initialized GuessSpecialSeqResult object 
          - ``combined_remseq`` -- The remaining sequence in the combined
                                   result 
          """
          AssertVariableType(match1, GuessSpecialSeqResult) 
          AssertVariableType(match2, GuessSpecialSeqResult) 
          spseq_factors = match1.get_spseq_factors() + match2.get_spseq_factors()
          combined_match = GuessSpecialSeqResult(combined_remseq, 
                                                 spseq_factors[0][0], 
                                                 spseq_factors[0][1])
          for (seqgen, lcidx_func) in spseq_factors[1:]: 
               combined_match.add_spseq_factor(seqgen, lcidx_func)
          ##
          return combined_match
     ##

     @property
     def remaining_sequence(self): 
          r""" 
          Property that returns the leftover, or remaining, sequence 
          after removing the matching special sequence factors. 
          """
          return self.remseq
     ## 

     @remaining_sequence.setter
     def remaining_sequence(self, remseq): 
          r""" 
          Sets the remaining_sequence property. 

          INPUTS: 
          - ``remseq`` -- A list of integer-like terms 
          """
          self.remseq = remseq
     ##

     def finalize_match(self): 
          r""" 
          Sets an internal marker indicating that the working match object 
          is complete and ready to be returned. (Convenience function for 
          building the complete matching formulas) 
          """ 
          self.finalize = True
     ##

     def is_final(self):
          r""" 
          Tests whether the match is finalized (i.e., whether finalize_match() 
          has been called on the object). 
          """ 
          return self.finalize
     ## 

     @staticmethod
     def exact_match(remseq): 
          r""" 
          Tests whether the list of special sequence factors in the formula 
          match exactly generate the input sequence. In other words, the 
          function tests whether the remaining sequence terms are all 
          identically equal to one. 

          INPUTS: 
          - ``remseq`` -- A remaining sequence (list of integer-like terms) 
          
          EXAMPLES: 

          sage: GuessSpecialSeqResult.exact_match([1, 1, 1, 1, 1, 1])
          True
          sage: GuessSpecialSeqResult.exact_match([1, 1, -1, 1, 2, 3])
          False
          
          """
          eqone_func = lambda e: e == 1
          remseq2 = filter(eqone_func, remseq)
          return len(remseq) == len(remseq2)
     ## 

     def remseq_has_match(self): 
          r""" 
          Tests whether we can obtain a comparitively simpler formula for the 
          remaining sequence terms by calling the existing helper methods 
          packaged in guess_remaining_sequence above. 
          """
          if GuessSpecialSeqResult.exact_match(self.remseq): 
               return True
          ##
          fmatches = guess_remaining_sequence(self.remseq)
          return len(fmatches) > 0
     ## 

     def get_spseq_factors(self): 
          r""" 
          Returns the local list of special sequence function factor data 
          in the form of tuples stored as 
          (SequenceGenerator, (LinearCombination, LinearCombination)). 
          """
          return self.spseq_factors
     ##

     def add_spseq_factor(self, seqgen, lcidx): 
          r""" 
          Adds a new special sequence factor to the list of special sequence 
          factors in the local match object. 

          INPUTS: 
          - ``seqgen`` -- A SequenceGenerator object for the special sequence 
          - ``lcidx`` -- A 2-tuple of LinearCombination objects representing
                         the inputs to the special sequence function when 
                         evaluated at the original sequence indices 
          
          See also get_spseq_factors above for the format these new data 
          items are stored in in the local instance. 

          """ 
          self.spseq_factors += [(seqgen, lcidx)]
     ## 

     def get_factor_string(self, seqgen, lcidx_func, domain_dim): 
          r""" 
          Helper method called by __str__ below which returns string 
          representations of the LinearCombination index objects for a 
          given set of input special sequence factor data. 

          INPUTS: 
          - ``seqgen`` -- SequenceGenerator object generating the given special
                          sequence factor 
          - ``lcidx_func`` -- 2-tuple of LinearCombination objects for the 
                              inputs to the given matching special sequence 
                              factor 
          - ``domain_dim`` -- The dimension of the domain of the given 
                              special sequence factor 
          """ 
          subst_args = []
          if domain_dim == 1: 
               subst_args = str(str(lcidx_func[X]) % 'n')
          else: 
               subst_args = (str(str(lcidx_func[X]) % 'n'), \
                             str(str(lcidx_func[Y]) % 'n'))
          ## 
          return str(seqgen) % subst_args
     ## 

     def __str__(self): 
          r""" 
          Returns a string representation of the matching sequence formula. 
          """
          rstr = ""
          for (ssidx, (seqgen, lcidx)) in enumerate(self.spseq_factors):
               if ssidx > 0: 
                    rstr += " * "
               ## 
               rstr += self.get_factor_string(seqgen, lcidx, \
                                              seqgen.domain_dimension())
          ## 
          rstr += " * RemSeq(n)"
          return rstr
     ##

     def print_summary(self): 
          r""" 
          Prints a summary of the local matching sequence formula. 
          """
          print "   >> Factors Formula:    ", str(self)
          print "   >> Remaining Sequence: ", self.remseq 
          n = var('n')
          remseq_formulas = guess_remaining_sequence(self.remseq)
          remseq_formulas = map(lambda f: f(n), remseq_formulas)
          print "   >> Remaining Sequence Formula(s): ", remseq_formulas
     ## 

## 

def guess_ssequence_factors(sequence, spec_seqs = SPECIAL_SEQUENCES, 
                            quick_eval = True): 
     r""" 
     Helper function that returns a list of formula matches corresponding to 
     special sequence factors of the input sequence with predictable 
     linear index combinations. 

     INPUTS: 
     - ``sequence`` -- An input sequence (list of integer-like values) 
     - ``spec_seqs`` -- A list of special SequenceGenerator objects 
     - ``quick_eval`` -- Whether to search for quick, rather than complete, 
                         matching formulas of the input sequence 
     """ 
     formula_matches = []
     for spseq in spec_seqs: 
          search_func = SearchSequenceFactors(seqgen = spseq)
          spfactors = search_func.compute_special_sequence_factors(sequence)
          for (lcidx, remseq) in spfactors: 
               fmatch = GuessSpecialSeqResult(remseq, spseq, lcidx)
               if quick_eval and fmatch.remseq_has_match():
                    fmatch.finalize_match()
               ##
               formula_matches += [fmatch]
          ##
     ## 
     return formula_matches
## 

def guess_special_sequence_factors(sequence, spseqs = SPECIAL_SEQUENCES, 
                                   user_spseqs = [], quick_eval = True): 
     r""" 

     INPUTS: 
     - ``sequence`` -- An input sequence (list of integer-like values)
     - ``spseqs`` -- A list of special SequenceGenerator objects
     - ``user_spseqs`` -- An auxiliary list of special SequenceGenerator objects
     - ``quick_eval`` -- Whether to search for quick, rather than complete, 
                         matching formulas of the input sequence 
     """ 
     init_guess_results = guess_remaining_sequence(sequence)
     special_seqs = user_spseqs + spseqs
     fmatches, working_fmatches = [], []
     working_fmatches += guess_ssequence_factors(sequence, 
                         special_seqs, quick_eval = quick_eval)
     
     temp_working_matches = []
     for (fmidx, fmatch) in enumerate(working_fmatches): 
          if fmatch.is_final(): 
               fmatches += [fmatch]
          else: 
               temp_working_matches += [fmatch]
          ##
     ##
     
     return init_guess_results + fmatches, temp_working_matches
## 

def guess(sequence, quick_eval = True, spseq_factors = SPECIAL_SEQUENCES, 
          num_seq_factors = 1, expected_factors = None, 
          user_guess_func = None, indeterminates = TODO(None)): 
     r""" 
     Wrapper method around the collection of formula guessing routines 
     defined above and in GuessingFormulasForSequences.py that takes 
     multiple configuration options to control the guessing parameters in 
     finding prescribed special sequence factors of the input sequence. 

     INPUTS: 
     - ``sequence`` -- A list of integer-like elements 
     - ``quick_eval`` -- Whether to stop the search process for additional 
                         special sequence factors once one complete formula 
                         exists for that particular match 
                         (of potentially many matches). 
                         This option should speed up the search procedure. 
     - ``spseq_factors`` -- A list of special SequenceGenerator objects 
     - ``num_seq_factors`` -- The maximum number of distinct special sequence 
                              factors to search for in the input sequence 
     - ``expected_factors`` -- A list of expected sequence factors in the 
                               input sequence list 
     - ``user_guess_func`` -- A single-input guess function providing 
                              user expected factors present in the sequence 
     - ``indeterminates`` -- Not currently implemented, reserved for future use  

     EXAMPLES: 

     We provide multiple examples to document the current functionality of the 
     guessing package routines. Let's first load packages and then 
     define the following helper function for the ANSI formatted pretty 
     printing of our formula results: 

     sage: from Guess import *
     sage: from SpecialSequences import *
     sage: from ColorizeString import print_color_string
     
     sage: def print_match_summary(fmatches): 
     ....:      for (fidx, fmatch) in enumerate(fmatches): 
     ....:           header_string = " ==== Formula Match # %d / %d ==== " % \
     ....:                           (fidx + 1, len(fmatches))
     ....:           print_color_string(header_string, 
     ....:                              opts = ['FGWHITE', 'BGBLACK', 'BOLD', 'UL'])
     ....:           print ""
     ....:           fmatch.print_summary()
     ....:      print "\n"
     
     To begin with, let's guess a formula for a 
     sequence that is already easily handled by existing software: 

     sage: sm = guess([0, 1, 4, 9, 16, 25, 36], spseq_factors = [])
     sage: print_match_summary(sm) 
      ==== Formula Match # 1 / 2 ==== 
         >> Sequence Formula:  n^2 + 2*n + 1
      ==== Formula Match # 2 / 2 ==== 
         >> Sequence Formula:  n^2 + 2*n + 1

     Another example: 
     
     sage: seq = [2**m for m in range(2, 10)]
     sage: print seq 
     [4, 8, 16, 32, 64, 128, 256, 512]
     sage: sm = guess(seq, spseq_factors = [])
     sage: sm[0].print_summary()
        >> Sequence Formula:  4*2^n

     Let's demonstrate another couple of options to the guess function 
     before moving on. First, we consider the usage of the ``expected_factors`` 
     option used to remove noisy special factors when we preprocess the 
     sequence: 
     
     sage: noise = [1, 11, 3, 17, 19, 18, 2, 10, 1]
     sage: seq = [(3**n) * noise[n] for n in range(0, len(noise))]
     sage: sm = guess(seq, spseq_factors = [], expected_factors = noise)
     sage: sm[0].print_summary()
        >> Sequence Formula:  3^n

     We can also define more predictable expected factors of the input 
     sequence by defining a user guess function: 

     sage: user_guess_func = lambda n: Rational(3 * (2**n)) / factorial(n)
     sage: seq = [3 * user_guess_func(n) * Binom2(n+3, 3) for n in range(1, 8)]
     sage: print seq 
     [9, 288, 1800, 4800, 7350, 37632/5, 28224/5, 23040/7]
     sage: sm = guess(seq, user_guess_func = user_guess_func, 
     ....:            spseq_factors = [SeqGen_Binomial, SeqGen_BinomialSquared])
     sage: print_match_summary(sm)
     TODO
     
     Next, lets define an input sequence involving the factorial normalized 
     first-order harmonic numbers. This sequence arises as cases of at least 
     three of the default special sequence factors, so this example 
     should generate some additional results for comparison: 

     sage: seq = [factorial(n+1) * harmonic_number(n+1, 1) for n in range(0, 8)]
     sage: print seq 
     sage: spseq_factors = [SeqGen_FirstOrderHarmonic, SeqGen_StirlingS1, 
     ....:                  SeqGen_NormalizedHarmonicNumber]
     sage: sm = guess(seq, spseq_factors = spseq_factors, num_seq_factors = 1)
     sage: print_summary_match(sm)
     TODO 

     Similarly, we repeat the same search except with index inputs to the 
     sequence shifted from n+1 -> 2n+1: 

     sage: seq = [factorial(2*n+1) * harmonic_number(2*n+1, 1) for n in range(0, 8)]
     sage: sm = guess(seq, spseq_factors = spseq_factors, num_seq_factors = 1)
     sage: print_summary_match(sm)
     TODO 

     Another related sequence case is formed by the second column of the 
     Stirling numbers of the first kind with S1(n, 2) = (n-1)! H_{n-1}, 
     where H_n denotes the sequence of first-order harmonic numbers 
     corresponding to the sequence of partial sums of the divergent 
     harmonic series: 


     MISC TODO / NOTES / FUTURE FEATURES: 
     - The FriCAS package returns an incorrect formula on the sequence 
       input of [1, -1, 1, -1, 1, -1, 1, -1]
     - In the tutorial, need to add more non-toy examples of special sequences 
       arising in applications 
     - Need to work on handling / optimizing cases of special sequence 
       factors with many entries of 1, such as in the Stirling number triangles. 
       This should greatly improve the orders of long running times for guess and 
       allow for plausible guessing over a larger number of special sequences at 
       once 
     - Would like to add polynomial sequence recognition to the package 
     - Add the GFun.py routines which include many sequence (generating
       function) transformations to try out on the input sequences 
     - For logistical purposes, we only handle special sequences with 
       domain dimensions of 1 or 2. These two cases should cover almost all 
       non-parameterized special sequence cases which are not polynomial in 
       one or more indeterminates, x. 

     """

     if len(sequence) < MIN_SEQUENCE_ELEMENTS: 
          raise ValueError("Insufficient number of sequence elements " + \
                           "must pass >= %d values for guessing" % \
                           MIN_SEQUENCE_ELEMENTS) 
          return None
     ##

     ## preprocess the sequence: 
     for (sidx, s) in enumerate(sequence): 
          if expected_factors != None and not IsZero(expected_factors[sidx]):
               s = s / expected_factors[sidx]
          if user_guess_func != None and not IsZero(user_guess_func(sidx)): 
               s = s / user_guess_func(sidx)
          ## 
          sequence[sidx] = s
     ## 
     
     ## strip any leading zeros for ease of processing: 
     while sequence[0] == 0: 
          sequence = sequence[1:]
     ## 
     print sequence

     guess_subroutine = lambda sequence: \
          guess_special_sequence_factors(sequence, quick_eval = quick_eval, 
                                         spseqs = spseq_factors, user_spseqs = [])
     merge_func = lambda prevfm, nextfm: \
          GuessSpecialSeqResult.merge(prevfm, nextfm, nextfm.remaining_sequence)
     
     fmatches, working_fmatches, next_fmatches = [], [], [] 
     fmatches, working_fmatches = guess_subroutine(sequence) 

     ## search for multiple sequence factors if specified: 
     for sfidx in range(1, num_seq_factors): 
          temp_working_fmatches = []
          for wfm in working_fmatches: 
               remseq_terms = wfm.remaining_sequence()
               next_fmatches, next_working_fmatches = guess_subrotine(remseq_terms)
               for nwfm in next_working_fmatches: 
                    temp_working_fmatches += [merge_func(wfm, nwfm)] 
               for fm in next_fmatches:
                    fmatches += [merge_func(wfm, fm)] 
               ## 
          ## 
          working_fmatches = temp_working_fmatches
     ## 
     
     for wfm in working_fmatches: 
          if wfm.remseq_has_match(): 
               fmatches += [wfm] 
          ##
     ##

     return fmatches

## 

