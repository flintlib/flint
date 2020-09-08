.. _qqbar:

**qqbar.h** -- algebraic numbers represented by minimal polynomials
===============================================================================

A :type:`qqbar_t` represents a real or complex algebraic number
(an element of `\overline{\mathbb{Q}}`) by its unique reduced
minimal polynomial in `\mathbb{Z}[x]` and an isolating complex interval.
The precision of isolating intervals is maintained automatically to
ensure that all operations on :type:`qqbar_t` instances are exact.

This representation is useful for working with
individual algebraic numbers of moderate degree (up to 100, say).
Arithmetic in this representation is expensive: an arithmetic operation
on numbers of degrees *m* and *n* involves computing and then factoring an
annihilating polynomial of degree *mn* and potentially also performing
numerical root-finding. For doing repeated arithmetic, it is generally
more efficient to work with the :type:`ca_t` type in a fixed
number field.
The :type:`qqbar_t` type is used internally by the :type:`ca_t` type
to represent the embedding of number fields in `\mathbb{R}` or
`\mathbb{C}` and to decide predicates for algebraic numbers.


Types and macros
-------------------------------------------------------------------------------

.. type:: qqbar_struct

.. type:: qqbar_t

    A *qqbar_struct* consists of an *fmpz_poly_struct* and an *acb_struct*.
    A *qqbar_t* is defined as an array of length one of type
    *qqbar_struct*, permitting a *qqbar_t* to be passed by
    reference.

.. type:: qqbar_ptr

    Alias for ``qqbar_struct *``, used for *qqbar* vectors.

.. type:: qqbar_srcptr

    Alias for ``const qqbar_struct *``, used for *qqbar* vectors
    when passed as readonly input to functions.

.. macro:: QQBAR_POLY(x)

    Macro returning a pointer to the minimal polynomial of *x* which can be used as an *fmpz_poly_t*.

.. macro:: QQBAR_COEFFS(x)

    Macro returning a pointer to the array of *fmpz* coefficients of the
    minimal polynomial of *x*.

.. macro:: QQBAR_ENCLOSURE(x)

    Macro returning a pointer to the enclosure of *x* which can be used as an *acb_t*.

Memory management
-------------------------------------------------------------------------------

.. function:: void qqbar_init(qqbar_t res)

    Initializes the variable *res* for use, and sets its value to zero.

.. function:: void qqbar_clear(qqbar_t res)

    Clears the variable *res*, freeing or recycling its allocated memory.

.. function:: qqbar_ptr qqbar_vec_init(slong len)

    Returns a pointer to an array of *len* initialized *qqbar_struct*:s.

.. function:: void qqbar_vec_clear(qqbar_ptr vec, slong len)

    Clears all *len* entries in the vector *vec* and frees the
    vector itself.

Assignment
-------------------------------------------------------------------------------

.. function:: void qqbar_swap(qqbar_t x, qqbar_t y)

    Swaps the values of *x* and *y* efficiently.

.. function:: void qqbar_set(qqbar_t res, const qqbar_t x)
              void qqbar_set_si(qqbar_t res, slong x)
              void qqbar_set_ui(qqbar_t res, ulong x)
              void qqbar_set_fmpz(qqbar_t res, const fmpz_t x)
              void qqbar_set_fmpq(qqbar_t res, const fmpq_t x)

    Sets *res* to the value *x*.

.. function:: void qqbar_set_re_im(qqbar_t res, const qqbar_t x, const qqbar_t y)

    Sets *res* to the value `x + yi`.

.. function:: int qqbar_set_d(qqbar_t res, double x)
              int qqbar_set_re_im_d(qqbar_t res, double x, double y)

    Sets *res* to the value *x* or `x + yi` respectively. These functions
    performs error handling: if *x* and *y* are finite, the conversion succeeds
    and the return flag is 1. If *x* or *y* is non-finite (infinity or NaN),
    the conversion fails and the return flag is 0.

Properties
-------------------------------------------------------------------------------

.. function:: slong qqbar_degree(const qqbar_t x)

    Returns the degree of *x*, i.e. the degree of the minimal polynomial.

.. function:: int qqbar_is_rational(const qqbar_t x)

    Returns whether *x* is a rational number.

.. function:: int qqbar_is_integer(const qqbar_t x)

    Returns whether *x* is an integer (an element of `\mathbb{Z}`).

.. function:: int qqbar_is_algebraic_integer(const qqbar_t x)

    Returns whether *x* is an algebraic integer, i.e. whether its minimal
    polynomial has leading coefficient 1.

.. function:: int qqbar_is_zero(const qqbar_t x)
              int qqbar_is_one(const qqbar_t x)
              int qqbar_is_neg_one(const qqbar_t x)

    Returns whether *x* is the number `0`, `1`, `-1`.

.. function:: int qqbar_is_i(const qqbar_t x)
              int qqbar_is_neg_i(const qqbar_t x)

    Returns whether *x* is the imaginary unit `i` (respectively `-i`).

.. function:: int qqbar_is_real(const qqbar_t x)

    Returns whether *x* is a real number.

.. function:: void qqbar_height(fmpz_t res, const qqbar_t x)

    Sets *res* to the height of *x* (the largest absolute value of the
    coefficients of the minimal polynomial of *x*).

.. function:: slong qqbar_height_bits(const qqbar_t x)

    Returns the height of *x* (the largest absolute value of the
    coefficients of the minimal polynomial of *x*) measured in bits.

.. function:: int qqbar_within_limits(const qqbar_t x, slong deg_limit, slong bits_limit)

    Checks if *x* has degree bounded by *deg_limit* and height
    bounded by *bits_limit* bits, returning 0 (false) or 1 (true).
    If *deg_limit* is set to 0, the degree check is skipped,
    and similarly for *bits_limit*.

.. function:: int qqbar_binop_within_limits(const qqbar_t x, const qqbar_t y, slong deg_limit, slong bits_limit)

    Checks if `x + y`, `x - y`, `x \cdot y` and `x / y` certainly have
    degree bounded by *deg_limit* (by multiplying the degrees for *x* and *y*
    to obtain a trivial bound). For *bits_limits*, the sum of the bit heights
    of *x* and *y* is checked against the bound (this is only a heuristic).
    If *deg_limit* is set to 0, the degree check is skipped,
    and similarly for *bits_limit*.


Special values
-------------------------------------------------------------------------------

.. function:: void qqbar_zero(qqbar_t res)

    Sets *res* to the number 0.

.. function:: void qqbar_one(qqbar_t res)

    Sets *res* to the number 1.

.. function:: void qqbar_i(qqbar_t res)

    Sets *res* to the imaginary unit `i`.

.. function:: void qqbar_phi(qqbar_t res)

    Sets *res* to the golden ratio `\varphi = \tfrac{1}{2}(\sqrt{5} + 1)`.

Input and output
-------------------------------------------------------------------------------

.. function:: void qqbar_print(const qqbar_t x)

    Prints *res* to standard output. The output shows the degree
    and the list of coefficients
    of the minimal polynomial followed by a decimal representation of
    the enclosing interval. This function is mainly intended for debugging.

.. function:: void qqbar_printn(const qqbar_t x, slong n)

    Prints *res* to standard output. The output shows a decimal
    approximation to *n* digits.

.. function:: void qqbar_printnd(const qqbar_t x, slong n)

    Prints *res* to standard output. The output shows a decimal
    approximation to *n* digits, followed by the degree of the number.

For example, *print*, *printn* and *printnd* with `n = 6` give
the following output for the numbers 0, 1, `i`, `\varphi`, `\sqrt{2}-\sqrt{3} i`:

.. code ::

    deg 1 [0, 1] 0
    deg 1 [-1, 1] 1.00000
    deg 2 [1, 0, 1] 1.00000*I
    deg 2 [-1, -1, 1] [1.61803398874989484820458683436563811772 +/- 6.00e-39]
    deg 4 [25, 0, 2, 0, 1] [1.4142135623730950488016887242096980786 +/- 8.67e-38] + [-1.732050807568877293527446341505872367 +/- 1.10e-37]*I

    0
    1.00000
    1.00000*I
    1.61803
    1.41421 - 1.73205*I

    0 (deg 1)
    1.00000 (deg 1)
    1.00000*I (deg 2)
    1.61803 (deg 2)
    1.41421 - 1.73205*I (deg 4)


Random generation
-------------------------------------------------------------------------------

.. function:: void qqbar_randtest(qqbar_t res, flint_rand_t state, slong deg, slong bits)

    Sets *res* to a random algebraic number with degree up to *deg* and
    with height (measured in bits) up to *bits*.

.. function:: void qqbar_randtest_real(qqbar_t res, flint_rand_t state, slong deg, slong bits)

    Sets *res* to a random real algebraic number with degree up to *deg* and
    with height (measured in bits) up to *bits*.

.. function:: void qqbar_randtest_nonreal(qqbar_t res, flint_rand_t state, slong deg, slong bits)

    Sets *res* to a random nonreal algebraic number with degree up to *deg* and
    with height (measured in bits) up to *bits*. Since all algebraic numbers
    of degree 1 are real, *deg* must be at least 2.

Comparisons
-------------------------------------------------------------------------------

.. function:: int qqbar_equal(const qqbar_t x, const qqbar_t y)

    Returns whether *x* and *y* are equal.

.. function:: int qqbar_cmp_re(const qqbar_t x, const qqbar_t y)

    Compares the real parts of *x* and *y*, returning -1, 0 or +1.

.. function:: int qqbar_cmp_im(const qqbar_t x, const qqbar_t y)

    Compares the imaginary parts of *x* and *y*, returning -1, 0 or +1.

.. function:: int qqbar_cmpabs_re(const qqbar_t x, const qqbar_t y)

    Compares the absolute values of the real parts of *x* and *y*, returning -1, 0 or +1.

.. function:: int qqbar_cmpabs_im(const qqbar_t x, const qqbar_t y)

    Compares the absolute values of the imaginary parts of *x* and *y*, returning -1, 0 or +1.

.. function:: int qqbar_cmpabs(const qqbar_t x, const qqbar_t y)

    Compares the absolute values of *x* and *y*, returning -1, 0 or +1.

Complex parts
-------------------------------------------------------------------------------

.. function:: void qqbar_conj(qqbar_t res, const qqbar_t x)

    Sets *res* to the complex conjugate of *x*.

.. function:: void qqbar_re(qqbar_t res, const qqbar_t x)

    Sets *res* to the real part of *x*.

.. function:: void qqbar_im(qqbar_t res, const qqbar_t x)

    Sets *res* to the imaginary part of *x*.

.. function:: void qqbar_re_im(qqbar_t res1, qqbar_t res2, const qqbar_t x)

    Sets *res1* to the real part of *x* and *res2* to the imaginary part of *x*.

.. function:: void qqbar_abs(qqbar_t res, const qqbar_t x)

    Sets *res* to the absolute value of *x*:

.. function:: void qqbar_abs2(qqbar_t res, const qqbar_t x)

    Sets *res* to the square of the absolute value of *x*.

.. function:: void qqbar_sgn(qqbar_t res, const qqbar_t x)

    Sets *res* to the complex sign of *x*, defined as 0 if *x* is zero
    and as `x / |x|` otherwise.

.. function:: int qqbar_sgn_re(const qqbar_t x)

    Returns the sign of the real part of *x* (-1, 0 or +1).

.. function:: int qqbar_sgn_im(const qqbar_t x)

    Returns the sign of the imaginary part of *x* (-1, 0 or +1).

.. function:: int qqbar_csgn(const qqbar_t x)

    Returns the extension of the real sign function taking the
    value 1 for *x* strictly in the right half plane, -1 for *x* strictly
    in the left half plane, and the sign of the imaginary part when *x* is on
    the imaginary axis. Equivalently, `\operatorname{csgn}(x) = x / \sqrt{x^2}`
    except that the value is 0 when *x* is zero.

Integer parts
-------------------------------------------------------------------------------

.. function:: void qqbar_floor(fmpz_t res, const qqbar_t x)

    Sets *res* to the floor function of *x*. If *x* is not real, the
    value is defined as the floor function of the real part of *x*.

.. function:: void qqbar_ceil(fmpz_t res, const qqbar_t x)

    Sets *res* to the ceiling function of *x*. If *x* is not real, the
    value is defined as the ceiling function of the real part of *x*.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void qqbar_neg(qqbar_t res, const qqbar_t x)

    Sets *res* to the negation of *x*.

.. function:: void qqbar_add(qqbar_t res, const qqbar_t x, const qqbar_t y)
              void qqbar_add_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
              void qqbar_add_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
              void qqbar_add_ui(qqbar_t res, const qqbar_t x, ulong y)
              void qqbar_add_si(qqbar_t res, const qqbar_t x, slong y)

    Sets *res* to the sum of *x* and *y*.

.. function:: void qqbar_sub(qqbar_t res, const qqbar_t x, const qqbar_t y)
              void qqbar_sub_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
              void qqbar_sub_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
              void qqbar_sub_ui(qqbar_t res, const qqbar_t x, ulong y)
              void qqbar_sub_si(qqbar_t res, const qqbar_t x, slong y)

    Sets *res* to the difference of *x* and *y*.

.. function:: void qqbar_mul(qqbar_t res, const qqbar_t x, const qqbar_t y)
              void qqbar_mul_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
              void qqbar_mul_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
              void qqbar_mul_ui(qqbar_t res, const qqbar_t x, ulong y)
              void qqbar_mul_si(qqbar_t res, const qqbar_t x, slong y)

    Sets *res* to the product of *x* and *y*.

.. function:: void qqbar_mul_2exp_si(qqbar_t res, const qqbar_t x, slong e)

    Sets *res* to *x* multiplied by `2^e`.

.. function:: void qqbar_sqr(qqbar_t res, const qqbar_t x)

    Sets *res* to the square of *x*.

.. function:: void qqbar_inv(qqbar_t res, const qqbar_t x, const qqbar_t y)

    Sets *res* to the multiplicative inverse of *y*.
    Division by zero calls *flint_abort*.

.. function:: void qqbar_div(qqbar_t res, const qqbar_t x, const qqbar_t y)
              void qqbar_div_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
              void qqbar_div_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
              void qqbar_div_ui(qqbar_t res, const qqbar_t x, ulong y)
              void qqbar_div_si(qqbar_t res, const qqbar_t x, slong y)
              void qqbar_fmpq_div(qqbar_t res, const fmpq_t x, const qqbar_t y)
              void qqbar_fmpz_div(qqbar_t res, const fmpz_t x, const qqbar_t y)
              void qqbar_ui_div(qqbar_t res, ulong x, const qqbar_t y)
              void qqbar_si_div(qqbar_t res, slong x, const qqbar_t y)

    Sets *res* to the quotient of *x* and *y*.
    Division by zero calls *flint_abort*.

.. function:: void qqbar_scalar_op(qqbar_t res, const qqbar_t x, const fmpz_t a, const fmpz_t b, const fmpz_t c)

    Sets *res* to the rational affine transformation `(ax+b)/c`, performed as
    a single operation. There are no restrictions on *a*, *b* and *c*
    except that *c* must be nonzero. Division by zero calls *flint_abort*.

Powers and roots
-------------------------------------------------------------------------------

.. function:: void qqbar_sqrt(qqbar_t res, const qqbar_t x)
              void qqbar_sqrt_ui(qqbar_t res, ulong x)

    Sets *res* to the principal square root of *x*.

.. function:: void qqbar_rsqrt(qqbar_t res, const qqbar_t x)

    Sets *res* to the reciprocal of the principal square root of *x*.
    Division by zero calls *flint_abort*.

.. function:: void qqbar_pow_ui(qqbar_t res, const qqbar_t x, ulong n)

    Sets *res* to *x* raised to the *n*-th power.

.. function:: void qqbar_root_ui(qqbar_t res, const qqbar_t x, ulong n)
              void qqbar_fmpq_root_ui(qqbar_t res, const fmpq_t x, ulong n)

    Sets *res* to the principal *n*-th root of *x*. The order *n*
    must be positive.

.. function:: void qqbar_fmpq_pow_si_ui(qqbar_t res, const fmpq_t x, slong m, ulong n)

    Sets *res* to the principal branch of `x^{m/n}`. The order *n*
    must be positive. Division by zero calls *flint_abort*.


Numerical enclosures
-------------------------------------------------------------------------------

The following functions guarantee a polished output in which both the real
and imaginary parts are accurate to *prec* bits and exact when exactly
representable (that is, when a real or imaginary part is a sufficiently
small dyadic number).
In some cases, the computations needed to polish the output may be
expensive. When polish is unnecessary, :func:`qqbar_enclosure_raw`
may be used instead. Alternatively, :func:`qqbar_cache_enclosure`
can be used to avoid recomputations.

.. function:: void qqbar_get_acb(acb_t res, const qqbar_t x, slong prec)

    Sets *res* to an enclosure of *x* rounded to *prec* bits.

.. function:: void qqbar_get_arb(arb_t res, const qqbar_t x, slong prec)

    Sets *res* to an enclosure of *x* rounded to *prec* bits, assuming that
    *x* is a real number. If *x* is not real, *res* is set to
    `[\operatorname{NaN} \pm \infty]`.

.. function:: void qqbar_get_arb_re(arb_t res, const qqbar_t x, slong prec)

    Sets *res* to an enclosure of the real part of *x* rounded to *prec* bits.

.. function:: void qqbar_get_arb_im(arb_t res, const qqbar_t x, slong prec)

    Sets *res* to an enclosure of the imaginary part of *x* rounded to *prec* bits.

.. function:: void qqbar_cache_enclosure(qqbar_t res, slong prec)

    Polishes the internal enclosure of *res* to at least *prec* bits
    of precision in-place. Normally, *qqbar* operations that need
    high-precision enclosures compute them on the fly without caching the results;
    if *res* will be used as an invariant operand for many operations,
    calling this function as a precomputation step can improve performance.


Conjugates
-------------------------------------------------------------------------------

.. function:: void qqbar_conjugates(qqbar_ptr res, const qqbar_t x)

    Sets the entries of the vector *res* to the *d* algebraic conjugates of
    *x*, including *x* itself, where *d* is the degree of *x*. The output is
    not guaranteed to be sorted in any particular order.

Polynomial evaluation
-------------------------------------------------------------------------------

.. function:: void _qqbar_evaluate_fmpq_poly(qqbar_t res, const fmpz * poly, const fmpz_t den, slong len, const qqbar_t x)
              void qqbar_evaluate_fmpq_poly(qqbar_t res, const fmpq_poly_t poly, const qqbar_t x)
              void _qqbar_evaluate_fmpz_poly(qqbar_t res, const fmpz * poly, slong len, const qqbar_t x)
              void qqbar_evaluate_fmpz_poly(qqbar_t res, const fmpz_poly_t poly, const qqbar_t x)

    Sets *res* to the value of the given polynomial *poly* evaluated at
    the algebraic number *x*. These methods detect simple special cases and
    automatically reduce *poly* if its degree is greater or equal
    to that of the minimal polynomial of *x*. In the generic case, evaluation
    simply uses Horner's rule.

.. function:: int qqbar_evaluate_fmpz_mpoly_iter(qqbar_t res, const fmpz_mpoly_t poly, qqbar_srcptr x, slong deg_limit, slong bits_limit, const fmpz_mpoly_ctx_t ctx)
              int qqbar_evaluate_fmpz_mpoly_horner(qqbar_t res, const fmpz_mpoly_t poly, qqbar_srcptr x, slong deg_limit, slong bits_limit, const fmpz_mpoly_ctx_t ctx)
              int qqbar_evaluate_fmpz_mpoly(qqbar_t res, const fmpz_mpoly_t poly, qqbar_srcptr x, slong deg_limit, slong bits_limit, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the value of *poly* evaluated at the algebraic numbers
    given in the vector *x*. The number of variables is defined by
    the context object *ctx*.

    The parameters *deg_limit* and *bits_limit*
    define evaluation limits: if any temporary result exceeds these limits
    (not necessarily the final value, in case of cancellation), the
    evaluation is aborted and 0 (failure) is returned. If evaluation
    succeeds, 1 is returned.

    The *iter* version iterates over all terms in succession and computes
    the powers that appear. The *horner* version uses a multivariate
    implementation of the Horner scheme. The default algorithm currently
    uses the Horner scheme.

Polynomial roots
-------------------------------------------------------------------------------

.. function:: void qqbar_roots_fmpz_poly(qqbar_ptr res, const fmpz_poly_t poly, int flags)
              void qqbar_roots_fmpq_poly(qqbar_ptr res, const fmpq_poly_t poly, int flags)

    Sets the entries of the vector *res* to the *d* roots of the polynomial
    *poly*. Roots with multiplicity appear with repetition in the
    output array.
    The output is not guaranteed to be sorted in any particular order,
    except that all instances of a repeated root always appear
    consecutively.

    The following *flags* are supported:

    - QQBAR_ROOTS_IRREDUCIBLE - if set, *poly* is assumed to be
      irreducible (it may still have constant content), and no polynomial
      factorization is performed internally.

.. function:: void qqbar_eigenvalues_fmpz_mat(qqbar_ptr res, const fmpz_mat_t mat, int flags)
              void qqbar_eigenvalues_fmpq_mat(qqbar_ptr res, const fmpz_mat_t mat, int flags)

    Sets the entries of the vector *res* to the eigenvalues of the
    square matrix *mat*. These functions compute the characteristic polynomial
    of *mat* and then call :func:`qqbar_roots_fmpz_poly` with the same
    flags.

Roots of unity and trigonometric functions
-------------------------------------------------------------------------------

The following functions use word-size integers *p* and *q*
instead of *fmpq_t* instances to express rational numbers.
This is to emphasize that
the computations are feasible only with small *q* in this representation
of algebraic numbers since the
associated minimal polynomials have degree `O(q)`.
The input *p* and *q* do not need to be reduced *a priori*,
but should not be close to the word boundaries (they may be added
and subtracted internally).

.. function:: void qqbar_root_of_unity(qqbar_t res, slong p, ulong q)

    Sets *res* to the root of unity `e^{2 \pi i p / q}`.

.. function:: int qqbar_is_root_of_unity(slong * p, ulong * q, const qqbar_t x)

    If *x* is not a root of unity, returns 0.
    If *x* is a root of unity, returns 1.
    If *p* and *q* are not *NULL* and *x* is a root of unity,
    this also sets *p* and *q* to the minimal integers with `0 \le p < q`
    such that `x = e^{2 \pi i p / q}`.

.. function:: void qqbar_exp_pi_i(qqbar_t res, slong p, ulong q)

    Sets *res* to the root of unity `e^{\pi i p / q}`.

.. function:: void qqbar_cos_pi(qqbar_t res, slong p, ulong q)
              void qqbar_sin_pi(qqbar_t res, slong p, ulong q)
              void qqbar_tan_pi(qqbar_t res, slong p, ulong q)
              void qqbar_cot_pi(qqbar_t res, slong p, ulong q)
              void qqbar_sec_pi(qqbar_t res, slong p, ulong q)
              void qqbar_csc_pi(qqbar_t res, slong p, ulong q)

    Sets *res* to the trigonometric function `\cos(\pi x)`,
    `\sin(\pi x)`, etc., with `x = \tfrac{p}{q}`.
    Evaluation at a pole of tan, cot, sec or csc inflicts division by zero.

.. function:: int qqbar_log_pi_i(slong * p, ulong * q, const qqbar_t x)

    If `y = \operatorname{log}(x) / (\pi i)` is algebraic, and hence
    necessarily rational, sets `y = p / q` to the reduced such
    fraction with `-1 < y \le 1` and returns 1.
    If *y* is not algebraic, returns 0.

.. function:: int qqbar_atan_pi(slong * p, ulong * q, const qqbar_t x)

    If `y = \operatorname{atan}(x) / \pi` is algebraic, and hence
    necessarily rational, sets `y = p / q` to the reduced such
    fraction with `|y| < \tfrac{1}{2}` and returns 1.
    If *y* is not algebraic, returns 0.

.. function:: int qqbar_asin_pi(slong * p, ulong * q, const qqbar_t x)

    If `y = \operatorname{asin}(x) / \pi` is algebraic, and hence
    necessarily rational, sets `y = p / q` to the reduced such
    fraction with `|y| \le \tfrac{1}{2}` and returns 1.
    If *y* is not algebraic, returns 0.

.. function:: int qqbar_acos_pi(slong * p, ulong * q, const qqbar_t x)

    If `y = \operatorname{acos}(x) / \pi` is algebraic, and hence
    necessarily rational, sets `y = p / q` to the reduced such
    fraction with `0 \le y \le 1` and returns 1.
    If *y* is not algebraic, returns 0.

.. function:: int qqbar_acot_pi(slong * p, ulong * q, const qqbar_t x)

    If `y = \operatorname{acot}(x) / \pi` is algebraic, and hence
    necessarily rational, sets `y = p / q` to the reduced such
    fraction with `-\tfrac{1}{2} < y \le \tfrac{1}{2}` and returns 1.
    If *y* is not algebraic, returns 0.

.. function:: int qqbar_asec_pi(slong * p, ulong * q, const qqbar_t x)

    If `y = \operatorname{asec}(x) / \pi` is algebraic, and hence
    necessarily rational, sets `y = p / q` to the reduced such
    fraction with `0 \le y \le 1` and returns 1.
    If *y* is not algebraic, returns 0.

.. function:: int qqbar_acsc_pi(slong * p, ulong * q, const qqbar_t x)

    If `y = \operatorname{acsc}(x) / \pi` is algebraic, and hence
    necessarily rational, sets `y = p / q` to the reduced such
    fraction with `-\tfrac{1}{2} \le y \le \tfrac{1}{2}` and returns 1.
    If *y* is not algebraic, returns 0.


Guessing and simplification
-------------------------------------------------------------------------------

.. function:: int qqbar_guess(qqbar_t res, const acb_t z, slong max_deg, slong max_bits, int flags, slong prec)

    Attempts to find an algebraic number *res* of degree at most *max_deg* and
    height at most *max_bits* bits matching the numerical enclosure *z*.
    The return flag indicates success.
    This is only a heuristic method, and the return flag neither implies a
    rigorous proof that *res* is the correct result, nor a rigorous proof
    that no suitable algebraic number with the given *max_deg* and *max_bits*
    exists. (Proof of nonexistence could in principle be computed,
    but this is not yet implemented.)

    The working precision *prec* should normally be the same as the precision
    used to compute *z*. It does not make much sense to run this algorithm
    with precision smaller than O(*max_deg* · *max_bits*).

    This function does a single iteration at the target *max_deg*, *max_bits*,
    and *prec*. For best performance, one should invoke this function
    repeatedly with successively larger parameters when the size of the
    intended solution is unknown or may be much smaller than a worst-case bound.

.. function:: int qqbar_express_in_field(fmpq_poly_t res, const qqbar_t alpha, const qqbar_t x, slong max_bits, int flags, slong prec)

    Attempts to express *x* in the number field generated by *alpha*, returning
    success (0 or 1). On success, *res* is set to a polynomial *f* of degree
    less than the degree of *alpha* and with height (counting both the numerator
    and the denominator, when the coefficients of *g* are
    put on a common denominator) bounded by *max_bits* bits, such that
    `f(\alpha) = x`.

    (Exception: the *max_bits* parameter is currently ignored if *x* is
    rational, in which case *res* is just set to the value of *x*.)

    This function looks for a linear relation heuristically using a working
    precision of *prec* bits. If *x* is expressible in terms of *alpha*,
    then this function is guaranteed to succeed when *prec* is taken
    large enough. The identity `f(\alpha) = x` is checked
    rigorously, i.e. a return value of 1 implies a proof of correctness.
    In principle, choosing a sufficiently large *prec* can be used to
    prove that *x* does not lie in the field generated by *alpha*,
    but the present implementation does not support doing so automatically.

    This function does a single iteration at the target *max_bits* and
    and *prec*. For best performance, one should invoke this function
    repeatedly with successively larger parameters when the size of the
    intended solution is unknown or may be much smaller than a worst-case bound.

Internal functions
-------------------------------------------------------------------------------

.. function:: void qqbar_fmpz_poly_composed_op(fmpz_poly_t res, const fmpz_poly_t A, const fmpz_poly_t B, int op)

    Given nonconstant polynomials *A* and *B*, sets *res* to a polynomial
    whose roots are `a+b`, `a-b`, `ab` or `a/b` for all roots *a* of *A*
    and all roots *b* of *B*. The parameter *op* selects the arithmetic
    operation: 0 for addition, 1 for subtraction, 2 for multiplication
    and 3 for division. If *op* is 3, *B* must not have zero as a root.

.. function:: void qqbar_binary_op(qqbar_t res, const qqbar_t x, const qqbar_t y, int op)

    Performs a binary operation using a generic algorithm. This does not
    check for special cases.

.. function:: int _qqbar_validate_uniqueness(acb_t res, const fmpz_poly_t poly, const acb_t z, slong max_prec)

    Given *z* known to be an enclosure of at least one root of *poly*,
    certifies that the enclosure contains a unique root, and in that
    case sets *res* to a new (possibly improved) enclosure for the same
    root, returning 1. Returns 0 if uniqueness cannot be certified.

    The enclosure is validated by performing a single step with the
    interval Newton method. The working precision is determined from the
    accuracy of *z*, but limited by *max_prec* bits.

    This method slightly inflates the enclosure *z* to improve the chances
    that the interval Newton step will succeed. Uniqueness on this larger
    interval implies uniqueness of the original interval, but not
    existence; when existence has not been ensured a priori,
    :func:`_qqbar_validate_existence_uniqueness` should be used instead.

.. function:: int _qqbar_validate_existence_uniqueness(acb_t res, const fmpz_poly_t poly, const acb_t z, slong max_prec)

    Given any complex interval *z*, certifies that the enclosure contains a
    unique root of *poly*, and in that case sets *res* to a new (possibly
    improved) enclosure for the same root, returning 1. Returns 0 if
    existence and uniqueness cannot be certified.

    The enclosure is validated by performing a single step with the
    interval Newton method. The working precision is determined from the
    accuracy of *z*, but limited by *max_prec* bits.

.. function:: void _qqbar_enclosure_raw(acb_t res, const fmpz_poly_t poly, const acb_t z, slong prec)
              void qqbar_enclosure_raw(acb_t res, const qqbar_t x, slong prec)

    Sets *res* to an enclosure of *x* accurate to about *prec* bits
    (the actual accuracy can be slightly lower, or higher).

    This function uses repeated interval Newton steps to polish the initial
    enclosure *z*, doubling the working precision each time. If any step
    fails to improve the accuracy significantly, the root is recomputed
    from scratch to higher precision.

    If the initial enclosure is accurate enough, *res* is set to this value
    without rounding and without further computation.

.. function:: int _qqbar_acb_lindep(fmpz * rel, acb_srcptr vec, slong len, int check, slong prec)

    Attempts to find an integer vector *rel* giving a linear relation between
    the elements of the real or complex vector *vec*, using the LLL algorithm.

    The working precision is set to the minimum of *prec* and the relative
    accuracy of *vec* (that is, the difference between the largest magnitude
    and the largest error magnitude within *vec*). 95% of the bits within the
    working precision are used for the LLL matrix, and the remaining 5% bits
    are used to validate the linear relation by evaluating the linear
    combination and checking that the resulting interval contains zero.
    This validation does not prove the existence or nonexistence
    of a linear relation, but it provides a quick heuristic way to eliminate
    spurious relations.

    If *check* is set, the return value indicates whether the validation
    was successful; otherwise, the return value simply indicates whether
    the algorithm was executed normally (failure may occur, for example,
    if the input vector is non-finite).

    In principle, this method can be used to produce a proof that no linear
    relation exists with coefficients up to a specified bit size, but this has
    not yet been implemented.


.. raw:: latex

    \newpage

