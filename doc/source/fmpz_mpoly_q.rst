.. _fmpz-mpoly-q:

**fmpz_mpoly_q.h** -- multivariate rational functions over Q
===============================================================================

An :type:`fmpz_mpoly_q_t` represents an element of 
`\mathbb{Q}(x_1,\ldots,x_n)` for fixed *n* as a pair of Flint 
multivariate polynomials (:type:`fmpz_mpoly_t`).
Instances are always kept in canonical form by ensuring that the GCD
of numerator and denominator is 1 and that the coefficient
of the leading term of the denominator is positive.

The user must create a multivariate polynomial context
(:type:`fmpz_mpoly_ctx_t`) specifying the number of variables *n* and
the monomial ordering.


Types and macros
-------------------------------------------------------------------------------

.. type:: fmpz_mpoly_q_struct

.. type:: fmpz_mpoly_q_t

    An *fmpz_mpoly_q_struct* consists of a pair of *fmpz_mpoly_struct*:s.
    An *fmpz_mpoly_q_t* is defined as an array of length one of type
    *fmpz_mpoly_q_struct*, permitting an *fmpz_mpoly_q_t* to be passed by
    reference.

.. macro:: fmpz_mpoly_q_numref(x)

    Macro returning a pointer to the numerator of *x* which can be used as an *fmpz_mpoly_t*.

.. macro:: fmpz_mpoly_q_denref(x)

    Macro returning a pointer to the denominator of *x* which can be used as an *fmpz_mpoly_t*.


Memory management
-------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_q_init(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)

    Initializes *res* for use, and sets its value to zero.

.. function:: void fmpz_mpoly_q_clear(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)

    Clears *res*, freeing or recycling its allocated memory.


Assignment
-------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_q_swap(fmpz_mpoly_q_t x, fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)

    Swaps the values of *x* and *y* efficiently.

.. function:: void fmpz_mpoly_q_set(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_set_fmpq(fmpz_mpoly_q_t res, const fmpq_t x, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_set_fmpz(fmpz_mpoly_q_t res, const fmpz_t x, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_set_si(fmpz_mpoly_q_t res, slong x, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the value *x*.


Canonicalisation
-------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_q_canonicalise(fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

    Puts the numerator and denominator of *x* in canonical form by removing
    common content and making the leading term of the denominator positive.

.. function:: int fmpz_mpoly_q_is_canonical(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

    Returns whether *x* is in canonical form.

    In addition to verifying that the numerator and denominator
    have no common content and that the leading term of the denominator 
    is positive, this function checks that the denominator is nonzero and that
    the numerator and denominator have correctly sorted terms
    (these properties should normally hold; verifying them
    provides an extra consistency check for test code).

Properties
-------------------------------------------------------------------------------

.. function:: int fmpz_mpoly_q_is_zero(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

    Returns whether *x* is the constant 0.

.. function:: int fmpz_mpoly_q_is_one(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

    Returns whether *x* is the constant 1.


Special values
-------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_q_zero(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the constant 0.

.. function:: void fmpz_mpoly_q_one(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the constant 1.

.. function:: void fmpz_mpoly_q_gen(fmpz_mpoly_q_t res, slong i, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the generator `x_{i+1}`.
    Requires `0 \le i < n` where *n* is the number of variables of *ctx*.


Input and output
-------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_q_print_pretty(const fmpz_mpoly_q_t f, const char ** x, fmpz_mpoly_ctx_t ctx)

    Prints *res* to standard output. If *x* is not *NULL*, the strings in
    *x* are used as the symbols for the variables.


Random generation
-------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_q_randtest(fmpz_mpoly_q_t res, flint_rand_t state, slong length, mp_limb_t coeff_bits, slong exp_bound, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to a random rational function where both numerator and denominator
    have up to *length* terms, coefficients up to size *coeff_bits*, and
    exponents strictly smaller than *exp_bound*.


Comparisons
-------------------------------------------------------------------------------

.. function:: int fmpz_mpoly_q_equal(const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)

    Returns whether *x* and *y* are equal.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_q_neg(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the negation of *x*.

.. function:: void fmpz_mpoly_q_add(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_add_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_add_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_add_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the sum of *x* and *y*.

.. function:: void fmpz_mpoly_q_sub(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_sub_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_sub_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_sub_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the difference of *x* and *y*.

.. function:: void fmpz_mpoly_q_mul(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_mul_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_mul_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the product of *x* and *y*.

.. function:: void fmpz_mpoly_q_div(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_div_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_div_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_div_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the quotient of *x* and *y*.
    Division by zero calls *flint_abort*.

.. function:: void fmpz_mpoly_q_inv(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the inverse of *x*. Division by zero
    calls *flint_abort*.

Content
-------------------------------------------------------------------------------

.. function:: void _fmpz_mpoly_q_content(fmpz_t num, fmpz_t den, const fmpz_mpoly_t xnum, const fmpz_mpoly_t xden, const fmpz_mpoly_ctx_t ctx)

.. function:: void fmpz_mpoly_q_content(fmpq_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the content of the coefficients of *x*.


.. raw:: latex

    \newpage

