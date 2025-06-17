.. _fmpz-mod-mpoly-q:

**fmpz_mod_mpoly_q.h** -- multivariate rational functions over Z/mZ
===============================================================================

An :type:`fmpz_mod_mpoly_q_t` represents an element of 
`\mathbb{F}_m(x_1,\ldots,x_n)` for some prime number *m* and fixed *n* as a pair of FLINT 
multivariate polynomials (:type:`fmpz_mod_mpoly_t`).
Instances are always kept in canonical form by ensuring that the GCD
of numerator and denominator is 1 and then normalizing denominator to a monic polynomial.

The user must create a multivariate polynomial context
(:type:`fmpz_mod_mpoly_ctx_t`) specifying the prime number *m* for the field `\mathbb{F}_m` 
the number of variables *n* and the monomial ordering. The user is responsible
for verifying that *m* is a prime number;
if *m* is composite, undefined behaviour may occur.


Types and macros
-------------------------------------------------------------------------------

.. type:: fmpz_mod_mpoly_q_struct

.. type:: fmpz_mod_mpoly_q_t

    An *fmpz_mod_mpoly_q_struct* consists of a pair of *fmpz_mod_mpoly_struct*:s.
    An *fmpz_mod_mpoly_q_t* is defined as an array of length one of type
    *fmpz_mod_mpoly_q_struct*, permitting an *fmpz_mod_mpoly_q_t* to be passed by
    reference.

.. macro:: fmpz_mod_mpoly_q_numref(x)

    Macro returning a pointer to the numerator of *x* which can be used as an *fmpz_mod_mpoly_t*.

.. macro:: fmpz_mod_mpoly_q_denref(x)

    Macro returning a pointer to the denominator of *x* which can be used as an *fmpz_mod_mpoly_t*.


Memory management
-------------------------------------------------------------------------------

.. function:: void fmpz_mod_mpoly_q_init(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)

    Initializes *res* for use, and sets its value to zero.

.. function:: void fmpz_mod_mpoly_q_clear(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)

    Clears *res*, freeing or recycling its allocated memory.


Assignment
-------------------------------------------------------------------------------

.. function:: void fmpz_mod_mpoly_q_swap(fmpz_mod_mpoly_q_t x, fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)

    Swaps the values of *x* and *y* efficiently.

.. function:: void fmpz_mod_mpoly_q_set(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_q_set_fmpq(fmpz_mod_mpoly_q_t res, const fmpq_t x, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_set_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_t x, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_set_si(fmpz_mod_mpoly_q_t res, slong x, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to the value *x*.
    The *fmpq* version returns 1 if the denominator of *x* is invertible,
    otherwise returns 0.


Canonicalisation
-------------------------------------------------------------------------------

.. function:: void fmpz_mod_mpoly_q_canonicalise(fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)

    Puts the numerator and denominator of *x* in canonical form by removing
    common content and making the denominator monic.

.. function:: int fmpz_mod_mpoly_q_is_canonical(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)

    Returns whether *x* is in canonical form.

    In addition to verifying that the numerator and denominator
    have no common content and that the denominator 
    is monic, this function checks that the denominator is nonzero and that
    the numerator and denominator have correctly sorted terms
    (these properties should normally hold; verifying them
    provides an extra consistency check for test code).

Properties
-------------------------------------------------------------------------------

.. function:: int fmpz_mod_mpoly_q_is_zero(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)

    Returns whether *x* is the constant 0.

.. function:: int fmpz_mod_mpoly_q_is_one(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)

    Returns whether *x* is the constant 1.

.. function:: void fmpz_mod_mpoly_q_used_vars(int * used, const fmpz_mod_mpoly_q_t f, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_used_vars_num(int * used, const fmpz_mod_mpoly_q_t f, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_used_vars_den(int * used, const fmpz_mod_mpoly_q_t f, const fmpz_mod_mpoly_ctx_t ctx)

    For each variable, sets the corresponding entry in *used* to the
    boolean flag indicating whether that variable appears in the
    rational function (respectively its numerator or denominator).

Special values
-------------------------------------------------------------------------------

.. function:: void fmpz_mod_mpoly_q_zero(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to the constant 0.

.. function:: void fmpz_mod_mpoly_q_one(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to the constant 1.

.. function:: void fmpz_mod_mpoly_q_gen(fmpz_mod_mpoly_q_t res, slong i, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to the generator `x_{i+1}`.
    Requires `0 \le i < n` where *n* is the number of variables of *ctx*.


Input and output
-------------------------------------------------------------------------------

The variable strings in *x* start with the variable of most significance at index `0`. If *x* is ``NULL``, the variables are named ``x1``, ``x2``, etc.

.. function:: void fmpz_mod_mpoly_q_print_pretty(const fmpz_mod_mpoly_q_t f, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)

    Prints *res* to standard output. If *x* is not *NULL*, the strings in
    *x* are used as the symbols for the variables.

.. function:: char * fmpz_mod_mpoly_q_get_str_pretty(const fmpz_mod_mpoly_q_t f, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)

    Return a string, which the user is responsible for cleaning up, representing *f*, given an array of variable strings *x*.

.. function:: int fmpz_mod_mpoly_q_set_str_pretty(fmpz_mod_mpoly_q_t res, const char * s, const char ** x, fmpz_mod_mpoly_ctx_t ctx)

    Set *res* to the fraction in the null-terminated string *str* given an array *x* of variable strings.
    If parsing *str* fails, *res* is set to zero, and `-1` is returned. Otherwise, `0` is returned.
    The operations ``+``, ``-``, ``*``, and ``/`` are permitted along with integers and the variables in *x*.
    The character ``^`` must be immediately followed by the (integer) exponent.
    If division by zero occurs, parsing fails.

Random generation
-------------------------------------------------------------------------------

.. function:: void fmpz_mod_mpoly_q_randtest(fmpz_mod_mpoly_q_t res, flint_rand_t state, slong length, slong exp_bound, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to a random rational function where both numerator and denominator
    have up to *length* terms and exponents strictly smaller than *exp_bound*.


Comparisons
-------------------------------------------------------------------------------

.. function:: int fmpz_mod_mpoly_q_equal(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)

    Returns whether *x* and *y* are equal.


Arithmetic
-------------------------------------------------------------------------------

The functions below which coerce from an *fmpq* or divide by an integer type
perform error handling by returning 1 if the denominator is invertible
and 0 otherwise.

.. function:: void fmpz_mod_mpoly_q_neg(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to the negation of *x*.

.. function:: void fmpz_mod_mpoly_q_add(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_q_add_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_add_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_add_fmpz_mod(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_add_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong y, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to the sum of *x* and *y*. 

.. function:: void fmpz_mod_mpoly_q_sub(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_q_sub_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_sub_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_sub_fmpz_mod(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_sub_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong y, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to the difference of *x* and *y*.

.. function:: void fmpz_mod_mpoly_q_mul(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_mul_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_mul_fmpz_mod(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
              void fmpz_mod_mpoly_q_mul_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong y, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to the product of *x* and *y*.

.. function:: void fmpz_mod_mpoly_q_div(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_q_div_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_q_div_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
              int fmpz_mod_mpoly_q_div_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong y, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to the quotient of *x* and *y*.

.. function:: void fmpz_mod_mpoly_q_inv(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)

    Sets *res* to the inverse of *x*. Division by zero
    calls *flint_abort*.



.. raw:: latex

    \newpage

