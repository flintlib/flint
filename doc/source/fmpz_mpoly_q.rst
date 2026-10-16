.. _fmpz-mpoly-q:

**fmpz_mpoly_q.h** -- multivariate rational functions over Q
===============================================================================

An :type:`fmpz_mpoly_q_t` represents an element of 
`\mathbb{Q}(x_1,\ldots,x_n)` for fixed *n* as a pair of FLINT 
multivariate polynomials (:type:`fmpz_mpoly_t`).
Instances are always kept in canonical form by ensuring that the GCD
of numerator and denominator is 1 and that the coefficient
of the leading term of the denominator is positive.
All functions assume their inputs are in canonical form and produce
canonical output; this is not checked. Use
:func:`fmpz_mpoly_q_canonicalise` on any directly constructed value.

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

.. function:: void fmpz_mpoly_q_used_vars(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_used_vars_num(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_used_vars_den(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)

    For each variable, sets the corresponding entry in *used* to the
    boolean flag indicating whether that variable appears in the
    rational function (respectively its numerator or denominator).

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

The variable strings in *x* start with the variable of most significance at index `0`. If *x* is ``NULL``, the variables are named ``x1``, ``x2``, etc.

.. function:: void fmpz_mpoly_q_print_pretty(const fmpz_mpoly_q_t f, const char ** x, const fmpz_mpoly_ctx_t ctx)

    Prints *res* to standard output. If *x* is not *NULL*, the strings in
    *x* are used as the symbols for the variables.

.. function:: char * fmpz_mpoly_q_get_str_pretty(const fmpz_mpoly_q_t f, const char ** x, const fmpz_mpoly_ctx_t ctx)

    Return a string, which the user is responsible for cleaning up, representing *f*, given an array of variable strings *x*.

.. function:: int fmpz_mpoly_q_set_str_pretty(fmpz_mpoly_q_t res, const char * s, const char ** x, fmpz_mpoly_ctx_t ctx)

    Set *res* to the fraction in the null-terminated string *str* given an array *x* of variable strings.
    If parsing *str* fails, *res* is set to zero, and `-1` is returned. Otherwise, `0` is returned.
    The operations ``+``, ``-``, ``*``, and ``/`` are permitted along with integers and the variables in *x*.
    The character ``^`` must be immediately followed by the (integer) exponent.
    If division by zero occurs, parsing fails.

Random generation
-------------------------------------------------------------------------------

.. function:: void fmpz_mpoly_q_randtest(fmpz_mpoly_q_t res, flint_rand_t state, slong length, ulong coeff_bits, slong exp_bound, const fmpz_mpoly_ctx_t ctx)

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
              void fmpz_mpoly_q_content(fmpq_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the content of the coefficients of *x*.


Rational functions with factored denominators
-------------------------------------------------------------------------------

An :type:`fmpz_mpoly_q_factored_t` represents an element of
`\mathbb{Q}(x_1,\ldots,x_n)` as a fraction `N / D` where the numerator
`N` is stored as an expanded :type:`fmpz_mpoly_t` and the denominator
`D = c \, p_1^{e_1} \cdots p_k^{e_k}` is stored in factored form as an
:type:`fmpz_mpoly_factor_t`.

This representation is useful in applications where
denominators are products of a modest number of simple factors that occur
repeatedly. Since the irreducible factors of the denominator are known,
addition, subtraction and multiplication only require trial divisions
by these factors instead of polynomial GCDs, and denominators never need to
be expanded. The tradeoff is that division (and inversion) requires
factoring the numerator of the divisor.

The canonical form is defined as follows:

* The constant *c* is a positive integer and the field *constant_den*
  of the denominator is 1.

* Each `p_i` is a nonconstant irreducible polynomial in
  `\mathbb{Z}[x_1,\ldots,x_n]` (hence primitive) with positive
  leading coefficient, and each exponent `e_i` is positive.

* The factors `p_i` are distinct and sorted in increasing order
  according to :func:`fmpz_mpoly_cmp`.

* The numerator is not divisible by any `p_i`, and the content of
  the numerator is coprime to *c*. Equivalently, the numerator and
  the expanded denominator have no common factor in `\mathbb{Z}[x_1,\ldots,x_n]`.

* The zero element is represented with `N = 0` and `D = 1`.

The generic ring wrapper :func:`gr_ctx_init_fmpz_mpoly_q_factored`
provides access to this type through the :ref:`gr <gr>` interface.

The canonical form is unique, so that equality can be checked
by structural comparison.
All functions assume that their inputs are in canonical form and produce
canonical output; this is not checked.
Use :func:`fmpz_mpoly_q_factored_canonicalise` on any directly
constructed value.

Unless otherwise stated, all functions allow aliasing between inputs
and outputs.

.. type:: fmpz_mpoly_q_factored_struct

.. type:: fmpz_mpoly_q_factored_t

    An *fmpz_mpoly_q_factored_struct* consists of an *fmpz_mpoly_struct*
    (the numerator) and an *fmpz_mpoly_factor_struct* (the denominator).
    An *fmpz_mpoly_q_factored_t* is defined as an array of length one of type
    *fmpz_mpoly_q_factored_struct*, permitting an *fmpz_mpoly_q_factored_t*
    to be passed by reference.

.. macro:: fmpz_mpoly_q_factored_numref(x)

    Macro returning a pointer to the numerator of *x* which can be used as an *fmpz_mpoly_t*.

.. macro:: fmpz_mpoly_q_factored_denref(x)

    Macro returning a pointer to the denominator of *x* which can be used as an *fmpz_mpoly_factor_t*.

.. function:: void fmpz_mpoly_q_factored_init(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)

    Initializes *res* for use, and sets its value to zero.

.. function:: void fmpz_mpoly_q_factored_clear(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)

    Clears *res*, freeing or recycling its allocated memory.

.. function:: void fmpz_mpoly_q_factored_swap(fmpz_mpoly_q_factored_t x, fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)

    Swaps the values of *x* and *y* efficiently.

.. function:: void fmpz_mpoly_q_factored_set(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_set_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_t x, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_set_fmpq(fmpz_mpoly_q_factored_t res, const fmpq_t x, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_set_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_t x, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_set_si(fmpz_mpoly_q_factored_t res, slong x, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the value *x*.

.. function:: void fmpz_mpoly_q_factored_zero(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_one(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the constant 0 or 1.

.. function:: void fmpz_mpoly_q_factored_gen(fmpz_mpoly_q_factored_t res, slong i, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the generator `x_{i+1}`.
    Requires `0 \le i < n` where *n* is the number of variables of *ctx*.

.. function:: void fmpz_mpoly_q_factored_set_fmpz_mpoly_q(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the value *x*. This requires factoring the
    denominator of *x*.

.. function:: void fmpz_mpoly_q_factored_get_fmpz_mpoly_q(fmpz_mpoly_q_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the value *x*, expanding the denominator.

.. function:: void fmpz_mpoly_q_factored_canonicalise(fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

    Puts *x* in canonical form. The input may have an arbitrary nonzero
    integer constant and *constant_den* in the denominator, and the bases of the
    denominator may be arbitrary nonzero polynomials (not necessarily
    irreducible, distinct or sorted) with arbitrary integer exponents;
    bases with negative exponents are moved to the numerator.
    Each base with positive exponent is factored into irreducibles,
    after which common factors are removed from the numerator and denominator.
    Throws if the denominator is zero.

.. function:: int fmpz_mpoly_q_factored_is_canonical(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

    Returns whether *x* is in canonical form. This function
    factors each base of the denominator to verify irreducibility
    and is intended for testing purposes.

.. function:: int fmpz_mpoly_q_factored_is_zero(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)
              int fmpz_mpoly_q_factored_is_one(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

    Returns whether *x* is the constant 0 or 1.

.. function:: int fmpz_mpoly_q_factored_equal(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)

    Returns whether *x* and *y* are equal.

.. function:: void fmpz_mpoly_q_factored_print_pretty(const fmpz_mpoly_q_factored_t x, const char ** vars, const fmpz_mpoly_ctx_t ctx)
              char * fmpz_mpoly_q_factored_get_str_pretty(const fmpz_mpoly_q_factored_t x, const char ** vars, const fmpz_mpoly_ctx_t ctx)

    Prints *x* to standard output or returns a string, which the user is
    responsible for freeing with *flint_free*. The denominator is printed in factored
    form, e.g. ``(x1+1)/(2*(x2)*(x1-x2)^3)``.
    If *vars* is *NULL*, the variables are named ``x1``, ``x2``, etc.

.. function:: int fmpz_mpoly_q_factored_set_str_pretty(fmpz_mpoly_q_factored_t res, const char * s, const char ** vars, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the rational function in the null-terminated string *s*
    given an array *vars* of variable strings (if *vars* is *NULL*, the
    variables are named ``x1``, ``x2``, etc.).
    The syntax is the same as for :func:`fmpz_mpoly_q_set_str_pretty`;
    in particular, this parses both the output of
    :func:`fmpz_mpoly_q_factored_get_str_pretty` and that of
    :func:`fmpz_mpoly_q_get_str_pretty`.
    The expression is evaluated using the arithmetic of this type,
    so each division by a non-constant expression factors that expression.
    If parsing fails (including if a division by zero occurs),
    *res* is set to zero and `-1` is returned. Otherwise, `0` is returned.

.. function:: void fmpz_mpoly_q_factored_randtest(fmpz_mpoly_q_factored_t res, flint_rand_t state, slong length, ulong coeff_bits, slong exp_bound, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to a random rational function where the numerator and
    the polynomials multiplied together to form the denominator have
    up to *length* terms, coefficients up to size *coeff_bits*, and
    exponents strictly smaller than *exp_bound*.

.. function:: void fmpz_mpoly_q_factored_neg(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the negation of *x*.

.. function:: void fmpz_mpoly_q_factored_add(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_sub(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the sum or difference of *x* and *y*.
    The common denominator is formed by taking the maximum exponent of each
    irreducible factor. Only factors occurring with the same exponent
    in both denominators can cancel against the new numerator; these
    are removed by trial division.

.. function:: void fmpz_mpoly_q_factored_mul(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the product of *x* and *y*.
    Factors of the denominator of *x* not occurring in the denominator of *y*
    are removed from the numerator of *y* by trial division, and vice versa.

.. function:: void fmpz_mpoly_q_factored_div(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the quotient of *x* and *y*. This requires factoring
    the numerator of *y*. Throws if *y* is zero.

.. function:: void fmpz_mpoly_q_factored_inv(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the inverse of *x*. This requires factoring
    the numerator of *x* and expanding the denominator of *x*.
    Throws if *x* is zero.

.. function:: void fmpz_mpoly_q_factored_pow_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t e, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_pow_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong e, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to *x* raised to the power *e*, which may be negative.
    A negative exponent requires factoring the numerator of *x*.
    We define `0^0 = 1`. Throws if *x* is zero and *e* is negative.

.. function:: void fmpz_mpoly_q_factored_add_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_add_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_add_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_add_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_sub_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_sub_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_sub_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_sub_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the sum or difference of *x* and *c*.

.. function:: void fmpz_mpoly_q_factored_mul_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_mul_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_mul_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_mul_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the product of *x* and *c*. Multiplication by an
    integer or rational number only changes the numerator and the constant
    of the denominator. Multiplication by a polynomial removes
    the factors of the denominator of *x* from *c* by trial division.

.. function:: void fmpz_mpoly_q_factored_div_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_div_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_div_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, const fmpz_mpoly_ctx_t ctx)
              void fmpz_mpoly_q_factored_div_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the quotient of *x* and *c*. Division by a
    polynomial requires factoring *c*. Throws if *c* is zero.


.. raw:: latex

    \newpage

