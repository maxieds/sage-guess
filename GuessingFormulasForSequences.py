
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
from sage.symbolic.expression import _eval_on_operands as eval_on_operands
from numpy import unique

class GuessFormulaResult(object): 
     r"""
     Base class that packages the data for a sequence formula 
     guess obtained from the FriCAS guess package
     """ 

     def __init__(self, sage_expr): 
          r"""
          Initializes the sequence function guess with input match data.

          INPUT:
          - ``sage_expr`` -- a sage expression representing the guess function
          """
          if isinstance(sage_expr, tuple): 
               sage_expr = sage_expr[0]
          ##
          self.sage_expr = sage_expr
          n = var('n')
          self.functor = lambda x: sage_expr.subs(n = x)
          #self.functor = None
     ## 

     def __str__(self): 
          r"""
          Returns a string representation of the sequence formula
          """
          n = var('n')
          return str(self.functor(n))
     ## 

     def __latex__(self): 
          r"""
          Returns a (LaTeX) string representation of the sequence formula
          """
          return str(self)
     ## 

     def __call__(self, x): 
          r"""
          operator() passes to and returns input from the sequence formula
          """
          return self.functor(x)
     ## 

     def get_functor(self): 
          r"""
          Returns a lambda function representing the sequence formula
          """
          return self.functor
     ## 

     def desc(self): 
          r"""
          Returns a description of the matching sequence formula. 
          Should be overridden by subclasses. 
          """
          return "<No Description>"
     ## 

     def print_summary(self): 
          r"""
          Prints a summary of the sequence formula match. 
          Should be overridden by subclasses.
          """
          n = var('n')
          print "   >> Sequence Formula: ", self.sage_expr #self.functor(n)
     ##

## GuessFormulaResult

class GuessingFormulasForSeqs(object): 
     r"""
     Container class for the python wrappers around the FriCAS guess package.
     - Implements: guessRat, guessExpRat, guessBinRat
     """

     @staticmethod
     def fricas2sage(fricas_return): 
          r"""
          Converts output returned by FriCAS into a Sage expression

          INPUT:

          - ``fricas_return`` -- output of calling a fricas.* function
          """
          if len(fricas_return) > 0:
               return factor(fricas_return[0].sage())
          else:
               return None
          
     ##

     @staticmethod
     def guessRationalFunction(sequence): 
          r"""
          Guesses a rational function formula for the input sequence.

          INPUT:

          - ``sequence`` -- a list of integers
          """
          fricas_output = fricas.guessRat(sequence)
          return GuessingFormulasForSeqs.fricas2sage(fricas_output)
     ##

     @staticmethod
     def guessExponentialRationalFunction(sequence): 
          r"""
          Guesses an exponential-rational function formula for the input sequence.

          INPUT:

          - ``sequence`` -- a list of integers
          """
          fricas_output = fricas.guessExpRat(sequence)
          return GuessingFormulasForSeqs.fricas2sage(fricas_output)
     ##

     @staticmethod
     def guessBinomialRationalFunction(sequence): 
          r"""
          Guesses a binomial-rational function formula for the input sequence.

          INPUT:

          - ``sequence`` -- a list of integers
          """
          fricas_output = fricas.guessBinRat(sequence)
          return GuessingFormulasForSeqs.fricas2sage(fricas_output)
     ##

     @staticmethod
     def guess(sequence): 
          r"""
          Guesses a formula for an input sequence using a multi-functional approach.

          INPUT:

          - ``sequence`` -- a list of integers 
          
          EXAMPLES: 

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

