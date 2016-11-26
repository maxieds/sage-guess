#### Guess.py : 
#### The main collaborative sequence formula guessing routines in Sage
#### (and some newly added routines to handle factors of special sequences)
#### Author:  Maxie D. Schmidt
#### Created: 2016.11.24 

from sage.all import *
from GuessingFormulasForSequences import *
from SpecialSequences import *
from SequenceElement import AssertVariableType
from SearchSequenceFactors import *
from ColorizeString import print_color_string
from GuessConfig import *

class OEISFormula(GuessFormulaResult): 
     
     def __init__(self, oeis_search_result): 
          self._seq_number = oeis_search_result.id()
          self._name = oeis_search_result.name()
          self.functor = sequence_func(list(oeis_search_result.first_terms()))
     ## 

     @property
     def sequence_number(self): 
          return self._seq_number
     ## 

     @property
     def name(self): 
          return self._name
     ## 

     def __str__(self): 
          return "%s(n)" % self.sequence_number
     ## 

     def __eq__(self, rhs): 
          AssertVariableType(rhs, OEISFormula)
          return self.sequence_number == rhs.sequence_number
     ## 

     def __ne__(self, rhs): 
          AssertVariableType(rhs, OEISFormula)
          return self.sequence_number != rhs.sequence_number
     ## 

     def desc(self): 
          return "%s: %s" % (self.sequence_number, self.name)
     ## 

     def print_summary(self): 
          print "   >>", self.desc()
          print "   >> First Sequence Values: ", \
                [self.functor(x) for x in range(0, 10)]
     ## 

## 

def guess_oeis_formula(sequence): 
     search_results = oeis(sequence, max_results = MAX_OEIS_RESULTS)
     search_results = map(OEISFormula, search_results)
     return search_results
## 

class GeneratingFunctionMatch(GuessFormulaResult): 

     def __init__(self, sequence, gfmatch): 
          self.gfmatch = gfmatch
          self.functor = sequence_func(sequence)
     ## 

     def __str__(self): 
          return "[x^n] " + str(self.gfmatch.ogf())
     ## 

     def __latex__(self): 
          return "[x^n] \\left(" + str(self.gfmatch.ogf()) + "\\right)"
     ## 
     
     def print_summary(self): 
          print "   >> Sequence OGF: ", self.gfmatch.ogf(), " = ", \
                self.gfmatch.series(4)
          print "   >> Recurrence: ", self.gfmatch.recurrence_repr()
     ## 


## GeneratingFunctionMatch

def guess_generating_function(sequence): 
     C, x = CFiniteSequences(QQ, 'x').objgen()
     smatch = C.guess(sequence) 
     if smatch == 0: 
          return []
     else: 
          return [GeneratingFunctionMatch(sequence, smatch)]
     ## 
## 

def guess_remaining_sequence(sequence): 
     guess_pkg_results = GuessingFormulasForSeqs.guess(sequence)
     #oeis_results = guess_oeis_formula(sequence) ## TODO
     #guess_gf_results = guess_generating_function(sequence) ## TODO
     return guess_pkg_results #+ oeis_results + guess_gf_results
     #return []
## 

class GuessSpecialSeqResult(GuessFormulaResult): 

     def __init__(self, remseq, seqgen, lcidx_func): 
          self.spseq_factors = [(seqgen, lcidx_func)]
          self.remseq = remseq
          self.finalize = False
     ## 

     @property
     def remaining_sequence(self): 
          return self.remseq
     ## 

     @remaining_sequence.setter
     def remaining_sequence(self, remseq): 
          self.remseq = remseq
     ##

     def finalize_match(self): 
          self.finalize = True
     ##

     def is_final(self):
          return self.finalize
     ## 

     @staticmethod
     def exact_match(remseq): 
          eqone_func = lambda e: e == 1
          remseq2 = filter(eqone_func, remseq)
          return len(remseq) == len(remseq2)
     ## 

     def remseq_has_match(self): 
          if GuessSpecialSeqResult.exact_match(self.remseq): 
               return True
          ##
          fmatches = guess_remaining_sequence(self.remseq)
          return len(fmatches) > 0
     ## 

     def get_spseq_factors(self): 
          return self.spseq_factors
     ##

     def add_spseq_factor(self, seqgen, lcidx): 
          self.spseq_factors += [(seqgen, lcidx)]
     ## 

     def get_factor_string(self, seqgen, lcidx_func, domain_dim): 
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

     def __latex__(self): 
          return str(TODO(self))
     ## 

     def print_summary(self): 
          print "   >> Factors Formula:    ", str(self)
          print "   >> Remaining Sequence: ", self.remseq 
     ## 

## 

def guess_ssequence_factors(sequence, spec_seqs = SPECIAL_SEQUENCES, 
                            quick_eval = True): 
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
                                   user_spseqs = [], quick_eval = True, 
                                   num_seq_factors = 1): 
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
     working_fmatches = temp_working_matches
     for fmatch in working_fmatches: 
          if fmatch.remseq_has_match(): 
               fmatches += [fmatch]
          ##
     ##
     return init_guess_results + fmatches
## 

def guess(sequence, quick_eval = True, spseq_factors = SPECIAL_SEQUENCES, 
          num_seq_factors = TODO(2), expected_factors = None, 
          user_guess_func = None, indeterminates = TODO(None)): 

     ## preprocess the sequence: 
     for (sidx, s) in enumerate(sequence): 
          if expected_factors != None:
               s = s / expected_factors[sidx]
          if user_guess_func != None: 
               s = s / user_guess_func(sidx)
          ## 
          sequence[sidx] = s
     ## 
     
     ## strip any leading zeros for ease of processing: 
     while sequence[0] == 0: 
          sequence = sequence[1:]
     ## 
     
     return guess_special_sequence_factors(sequence, quick_eval = quick_eval, 
            spseqs = spseq_factors, user_spseqs = [], 
            num_seq_factors = num_seq_factors)
## 

