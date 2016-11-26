
r"""
Python wrappers around the FriCAS GuessingFormulasForSequences package.

AUTHORS:

- Maxie D. Schmidt

EXAMPLE::

    sage: from GuessingFormulasForSequences import *
    sage: sm = GuessingFormulasForSeqs.guess([1, 4, 9, 16, 25, 36, 49])
    sage: sm[0].print_summary()
       >> Sequence Formula:  n^2 + 2*n + 1
    sage: func = sm[0].get_functor()
    sage: n = var('n'); factor(func(n))
    (n + 1)^2
    sage: func(3)
    16

"""

from sage.all import *
from sage.symbolic.function_factory import eval_on_operands
from numpy import unique

class GuessFormulaResult(object): 
     r"""
     Abstract class that packages the data for a sequence formula 
     guess obtained from the FriCAS guess package
     """ 

     def __init__(self, sage_expr): 
          self.sage_expr = sage_expr
          self.functor = lambda x: \
               eval_on_operands(lambda x: sage_expr)(n).subs(n == x)
     ## 

     def __str__(self): 
          n = var('n')
          return str(self.functor(n))
     ## 

     def __latex__(self): 
          return str(self)
     ## 

     def __call__(self, x): 
          return self.functor(x)
     ## 

     def get_functor(self): 
          return self.functor
     ## 

     def desc(self): 
          return "<No Description>"
     ## 

     def print_summary(self): 
          n = var('n')
          print "   >> Sequence Formula: ", self.functor(n)
     ##

## GuessFormulaResult

class GuessingFormulasForSeqs(object): 
     r"""
     Container class for the python wrappers around the FriCAS guess package.
     - Implements: guessRat, guessExpRat, guessBinRat
     """

     @staticmethod
     def fricas2sage(fricas_return): 
          n = var('n')
          fricas_return = fricas_return.unparsed_input_form() 
          if fricas_return != '[]': 
               return sage_eval(fricas_return[1:-1], locals = {'n': n}) 
          else: 
               return None
     ##

     @staticmethod
     def guessRationalFunction(sequence): 
          fricas_output = fricas.guessRat(sequence)
          return GuessingFormulasForSeqs.fricas2sage(fricas_output)
     ##

     @staticmethod
     def guessExponentialRationalFunction(sequence): 
          fricas_output = fricas.guessExpRat(sequence)
          return GuessingFormulasForSeqs.fricas2sage(fricas_output)
     ##

     @staticmethod
     def guessBinomialRationalFunction(sequence): 
          fricas_output = fricas.guessBinRat(sequence)
          return GuessingFormulasForSeqs.fricas2sage(fricas_output)
     ##

     @staticmethod
     def guess(sequence): 
          guess_funcs = [ 
               GuessingFormulasForSeqs.guessRationalFunction, 
               GuessingFormulasForSeqs.guessExponentialRationalFunction, 
               GuessingFormulasForSeqs.guessBinomialRationalFunction, 
          ]
          guess_formulas = map(lambda func: func(sequence), guess_funcs)
          filter_func = lambda elem: elem != None
          guess_formulas = filter(filter_func, guess_formulas)
          guess_formulas = unique(guess_formulas)
          guess_formulas = map(GuessFormulaResult, guess_formulas)
          return guess_formulas
    ## 

## GuessingFormulasForSequences

