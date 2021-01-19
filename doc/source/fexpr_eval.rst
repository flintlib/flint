.. _fexpr-eval:

**fexpr_eval.h** -- evaluation of symbolic expressions
===============================================================================

Evaluation is the mapping `E \to R` where *E* is the set of symbolic
expressions and *R* is a concrete set of values. This module will
support evaluating expressions for various such sets *R*, using
Flint, Arb, Antic and Calcium types, including:

* Booleans
* Rational numbers
* Algebraic numbers
* Complex numbers
* Singularities
* Finite fields
* Multivariate polynomials and rational functions
* Matrices

Ideally, evaluation will consist of traversing an expression and
performing calculations using a concrete type realizing the set *R*.
Making things more difficult, evaluation may require passing
through several different types internally. It may also require
symbolic rewriting. This module is meant to provide routines
for handling such issues automatically.

Types and macros
-------------------------------------------------------------------------------

.. type:: fexpr_eval_struct

.. type:: fexpr_eval_t

    An *fexpr_eval_struct* defines a context for interpreting a
    symbolic expression. It may contain the following data:

    * Symbol bindings
    * Assumptions and inferences
    * Cached data
    * Evaluation options

    An *fexpr_eval_t* is defined as an array of length one of type
    *fexpr_eval_struct*, permitting an *fexpr_eval_t* to be passed by
    reference.

Evaluation functions return one of the following status codes.

.. macro:: FEXPR_EVAL_SUCCESS

    The evaluation completed successfully.

.. macro:: FEXPR_EVAL_ERR_NOT_IMPLEMENTED

    The evaluation encountered a case for which an algorithm has not
    been implemented.

.. macro:: FEXPR_EVAL_ERR_SYNTAX

    The expression contains a "syntax" error (e.g. assignment to
    a non-symbol).

.. macro:: FEXPR_EVAL_ERR_SYMBOL

    The expression contains an undefined symbol in a place where
    a constant value is expected.

.. macro:: FEXPR_EVAL_ERR_TYPE

    The expression contains a type error (e.g. a number passed as
    input to a function that expects a boolean).

.. macro:: FEXPR_EVAL_ERR_OVERFLOW

    Evaluation could not complete due to requiring the construction
    of a too large exact object.

.. macro:: FEXPR_EVAL_ERR_PRECISION

    Evaluation could not complete due to numerical precision being
    exhausted.

Evaluation context
-------------------------------------------------------------------------------

.. function:: void fexpr_eval_init(fexpr_eval_t eval)

    Initializes the evaluation context *eval* for use.

.. function:: void fexpr_eval_clear(fexpr_eval_t eval)

    Clears the evaluation context *eval*.

.. function:: void fexpr_eval_push_def(fexpr_eval_t eval, const fexpr_t symbol, const fexpr_t value)

    Pushes the definition `x = v` where *x* is given by *symbol*
    and *v* is given by *value*. As long as this definition is active,
    instances of *x* will be replaced by *v* during evaluation.

.. function:: void fexpr_eval_pop_def(fexpr_eval_t eval)

    Pops the last pushed definition.

.. function:: void fexpr_eval_push_assumptions(fexpr_eval_t eval, const fexpr_t assumptions)

    Pushes assumptions: we assert that *assumptions* is a true
    statement, which will be exploited during evaluation.
    Several assumptions can be pushed simultaneously
    with an ``And``-expression. Examples of assumptions:

    * ``Element(n, ZZ)`` (`n \in \mathbb{Z}`)
    * ``And(Element(x, RR), Greater(x, 0))`` (`x \in \mathbb{R} \operatorname{and} x > 0`)
    * ``NotEqual(x, y)`` (`x \ne y`)
    * ``RiemannHypothesis`` (the Riemann hypothesis is assumed to be true)

.. function:: void fexpr_eval_pop_assumptions(fexpr_eval_t eval)

    Pops the last pushed assumptions.


.. raw:: latex

    \newpage
