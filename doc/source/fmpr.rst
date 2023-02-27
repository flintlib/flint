.. _fmpr:

**fmpr.h** -- Arb 1.x floating-point numbers (deprecated)
===============================================================================

This module is deprecated, and any methods contained herein could disappear
in the future. This module is mainly kept for testing the faster
:type:`arf_t` type that was introduced in Arb 2.0. Please use
:type:`arf_t` instead of the :type:`fmpr_t` type.

A variable of type :type:`fmpr_t` holds an arbitrary-precision binary
floating-point number, i.e. a rational number of the form
`x \times 2^y` where `x, y \in \mathbb{Z}` and `x` is odd;
or one of the special values zero, plus infinity, minus infinity,
or NaN (not-a-number).

The component `x` is called the *mantissa*, and `y` is called the
*exponent*. Note that this is just one among many possible
conventions: the mantissa (alternatively *significand*) is
sometimes viewed as a fraction in the interval `[1/2, 1)`, with the
exponent pointing to the position above the top bit rather than the
position of the bottom bit, and with a separate sign.

The conventions for special values largely follow those of the
IEEE floating-point standard. At the moment, there is no support
for negative zero, unsigned infinity, or a NaN with a payload, though
some these might be added in the future.

An *fmpr* number is exact and has no inherent "accuracy". We
use the term *precision* to denote either the target precision of
an operation, or the bit size of a mantissa (which in general is
unrelated to the "accuracy" of the number: for example, the
floating-point value 1 has a precision of 1 bit in this sense and is
simultaneously an infinitely accurate approximation of the
integer 1 and a 2-bit accurate approximation of
`\sqrt 2 = 1.011010100\ldots_2`).

Except where otherwise noted, the output of an operation is the
floating-point number obtained by taking the inputs as exact numbers,
in principle carrying out the operation exactly, and rounding the
resulting real number to the nearest representable floating-point
number whose mantissa has at most the specified number of bits, in
the specified direction of rounding. Some operations are always
or optionally done exactly.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpr_struct

    An *fmpr_struct* holds a mantissa and an exponent.
    If the mantissa and exponent are sufficiently small, their values are
    stored as immediate values in the *fmpr_struct*; large values are
    represented by pointers to heap-allocated arbitrary-precision integers.
    Currently, both the mantissa and exponent are implemented using
    the FLINT *fmpz* type. Special values are currently encoded
    by the mantissa being set to zero.

.. type:: fmpr_t

    An *fmpr_t* is defined as an array of length one of type
    *fmpr_struct*, permitting an *fmpr_t* to be passed by
    reference.

.. type:: fmpr_rnd_t

    Specifies the rounding mode for the result of an approximate operation.

.. macro:: FMPR_RND_DOWN

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction towards zero.

.. macro:: FMPR_RND_UP

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction away from zero.

.. macro:: FMPR_RND_FLOOR

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction towards minus infinity.

.. macro:: FMPR_RND_CEIL

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction towards plus infinity.

.. macro:: FMPR_RND_NEAR

    Specifies that the result of an operation should be rounded to the
    nearest representable number, rounding to an odd mantissa if there is a tie
    between two values. *Warning*: this rounding mode is currently
    not implemented (except for a few conversions functions where this 
    stated explicitly).

.. macro:: FMPR_PREC_EXACT

    If passed as the precision parameter to a function, indicates that no
    rounding is to be performed. This must only be used when it is known
    that the result of the operation can be represented exactly and fits
    in memory (the typical use case is working small integer values).
    Note that, for example, adding two numbers whose exponents are far
    apart can easily produce an exact result that is far too large to
    store in memory.

Memory management
-------------------------------------------------------------------------------

.. function:: void fmpr_init(fmpr_t x)

    Initializes the variable *x* for use. Its value is set to zero.

.. function:: void fmpr_clear(fmpr_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.


Special values
-------------------------------------------------------------------------------

.. function:: void fmpr_zero(fmpr_t x)

.. function:: void fmpr_one(fmpr_t x)

.. function:: void fmpr_pos_inf(fmpr_t x)

.. function:: void fmpr_neg_inf(fmpr_t x)

.. function:: void fmpr_nan(fmpr_t x)

    Sets *x* respectively to 0, 1, `+\infty`, `-\infty`, NaN.

.. function:: int fmpr_is_zero(const fmpr_t x)

.. function:: int fmpr_is_one(const fmpr_t x)

.. function:: int fmpr_is_pos_inf(const fmpr_t x)

.. function:: int fmpr_is_neg_inf(const fmpr_t x)

.. function:: int fmpr_is_nan(const fmpr_t x)

    Returns nonzero iff *x* respectively equals
    0, 1, `+\infty`, `-\infty`, NaN.

.. function:: int fmpr_is_inf(const fmpr_t x)

    Returns nonzero iff *x* equals either `+\infty` or `-\infty`.

.. function:: int fmpr_is_normal(const fmpr_t x)

    Returns nonzero iff *x* is a finite, nonzero floating-point value, i.e.
    not one of the special values 0, `+\infty`, `-\infty`, NaN.

.. function:: int fmpr_is_special(const fmpr_t x)

    Returns nonzero iff *x* is one of the special values
    0, `+\infty`, `-\infty`, NaN, i.e. not a finite, nonzero
    floating-point value.

.. function:: int fmpr_is_finite(fmpr_t x)

    Returns nonzero iff *x* is a finite floating-point value,
    i.e. not one of the values `+\infty`, `-\infty`, NaN.
    (Note that this is not equivalent to the negation of
    :func:`fmpr_is_inf`.)


Assignment, rounding and conversions
-------------------------------------------------------------------------------

.. function:: slong _fmpr_normalise(fmpz_t man, fmpz_t exp, slong prec, fmpr_rnd_t rnd)

    Rounds the mantissa and exponent in-place.

.. function:: void fmpr_set(fmpr_t y, const fmpr_t x)

    Sets *y* to a copy of *x*.

.. function:: void fmpr_swap(fmpr_t x, fmpr_t y)

    Swaps *x* and *y* efficiently.

.. function:: slong fmpr_set_round(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_set_round_fmpz(fmpr_t y, const fmpz_t x, slong prec, fmpr_rnd_t rnd)

    Sets *y* to a copy of *x* rounded in the direction specified by rnd to the
    number of bits specified by prec.

.. function:: slong _fmpr_set_round_mpn(slong * shift, fmpz_t man, mp_srcptr x, mp_size_t xn, int negative, slong prec, fmpr_rnd_t rnd)

    Given an integer represented by a pointer *x* to a raw array of
    *xn* limbs (negated if *negative* is nonzero), sets *man* to
    the corresponding floating-point mantissa rounded to *prec* bits in
    direction *rnd*, sets *shift* to the exponent, and returns the error bound.
    We require that *xn* is positive and that the leading limb of *x* is nonzero.

.. function:: slong fmpr_set_round_ui_2exp_fmpz(fmpr_t z, mp_limb_t lo, const fmpz_t exp, int negative, slong prec, fmpr_rnd_t rnd)

    Sets *z* to the unsigned integer *lo* times two to the power *exp*,
    negating the value if *negative* is nonzero, and rounding the result
    to *prec* bits in direction *rnd*.

.. function:: slong fmpr_set_round_uiui_2exp_fmpz(fmpr_t z, mp_limb_t hi, mp_limb_t lo, const fmpz_t exp, int negative, slong prec, fmpr_rnd_t rnd)

    Sets *z* to the unsigned two-limb integer *{hi, lo}* times two to the
    power *exp*, negating the value if *negative* is nonzero, and rounding the result
    to *prec* bits in direction *rnd*.

.. function:: void fmpr_set_error_result(fmpr_t err, const fmpr_t result, slong rret)

    Given the return value *rret* and output variable *result* from a
    function performing a rounding (e.g. *fmpr_set_round* or *fmpr_add*), sets
    *err* to a bound for the absolute error.

.. function:: void fmpr_add_error_result(fmpr_t err, const fmpr_t err_in, const fmpr_t result, slong rret, slong prec, fmpr_rnd_t rnd)

    Like *fmpr_set_error_result*, but adds *err_in* to the error.

.. function:: void fmpr_ulp(fmpr_t u, const fmpr_t x, slong prec)

    Sets *u* to the floating-point unit in the last place (ulp) of *x*.
    The ulp is defined as in the MPFR documentation and satisfies
    `2^{-n} |x| < u \le 2^{-n+1} |x|` for any finite nonzero *x*.
    If *x* is a special value, *u* is set to the absolute value of *x*.

.. function:: int fmpr_check_ulp(const fmpr_t x, slong r, slong prec)

    Assume that *r* is the return code and *x* is the floating-point result
    from a single floating-point rounding. Then this function returns nonzero
    iff *x* and *r* define an error of exactly 0 or 1 ulp.
    In other words, this function checks that :func:`fmpr_set_error_result`
    gives exactly 0 or 1 ulp as expected.

.. function:: int fmpr_get_mpfr(mpfr_t x, const fmpr_t y, mpfr_rnd_t rnd)

    Sets the MPFR variable *x* to the value of *y*. If the
    precision of *x* is too small to allow *y* to be represented
    exactly, it is rounded in the specified MPFR rounding mode.
    The return value indicates the direction of rounding,
    following the standard convention of the MPFR library.

.. function:: void fmpr_set_mpfr(fmpr_t x, const mpfr_t y)

    Sets *x* to the exact value of the MPFR variable *y*.

.. function:: double fmpr_get_d(const fmpr_t x, fmpr_rnd_t rnd)

    Returns *x* rounded to a *double* in the direction specified by *rnd*.

.. function:: void fmpr_set_d(fmpr_t x, double v)

    Sets *x* the the exact value of the argument *v* of type *double*.

.. function:: void fmpr_set_ui(fmpr_t x, ulong c)

.. function:: void fmpr_set_si(fmpr_t x, slong c)

.. function:: void fmpr_set_fmpz(fmpr_t x, const fmpz_t c)

    Sets *x* exactly to the integer *c*.

.. function:: void fmpr_get_fmpz(fmpz_t z, const fmpr_t x, fmpr_rnd_t rnd)

    Sets *z* to *x* rounded to the nearest integer in the direction
    specified by *rnd*. If *rnd* is *FMPR_RND_NEAR*, rounds to the
    nearest even integer in case of a tie.
    Aborts if *x* is infinite, NaN or if the exponent is unreasonably large.

.. function:: slong fmpr_get_si(const fmpr_t x, fmpr_rnd_t rnd)

    Returns *x* rounded to the nearest integer in the direction
    specified by *rnd*. If *rnd* is *FMPR_RND_NEAR*, rounds to the
    nearest even integer in case of a tie.
    Aborts if *x* is infinite, NaN, or the
    value is too large to fit in an *slong*.

.. function:: void fmpr_get_fmpq(fmpq_t y, const fmpr_t x)

    Sets *y* to the exact value of *x*. The result is undefined
    if *x* is not a finite fraction.

.. function:: slong fmpr_set_fmpq(fmpr_t x, const fmpq_t y, slong prec, fmpr_rnd_t rnd)

    Sets *x* to the value of *y*, rounded according to *prec* and *rnd*.

.. function:: void fmpr_set_fmpz_2exp(fmpr_t x, const fmpz_t man, const fmpz_t exp)

.. function:: void fmpr_set_si_2exp_si(fmpr_t x, slong man, slong exp)

.. function:: void fmpr_set_ui_2exp_si(fmpr_t x, ulong man, slong exp)

    Sets *x* to `\mathrm{man} \times 2^{\mathrm{exp}}`.

.. function:: slong fmpr_set_round_fmpz_2exp(fmpr_t x, const fmpz_t man, const fmpz_t exp, slong prec, fmpr_rnd_t rnd)

    Sets *x* to `\mathrm{man} \times 2^{\mathrm{exp}}`, rounded according
    to *prec* and *rnd*.

.. function:: void fmpr_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const fmpr_t x)

    Sets *man* and *exp* to the unique integers such that
    `x = \mathrm{man} \times 2^{\mathrm{exp}}` and *man* is odd,
    provided that *x* is a nonzero finite fraction.
    If *x* is zero, both *man* and *exp* are set to zero. If *x* is
    infinite or NaN, the result is undefined.

.. function:: int fmpr_get_fmpz_fixed_fmpz(fmpz_t y, const fmpr_t x, const fmpz_t e)

.. function:: int fmpr_get_fmpz_fixed_si(fmpz_t y, const fmpr_t x, slong e)

    Converts *x* to a mantissa with predetermined exponent, i.e. computes
    an integer *y* such that `y \times 2^e \approx x`, truncating if necessary.
    Returns 0 if exact and 1 if truncation occurred.


Comparisons
-------------------------------------------------------------------------------

.. function:: int fmpr_equal(const fmpr_t x, const fmpr_t y)

    Returns nonzero iff *x* and *y* are exactly equal. This function does
    not treat NaN specially, i.e. NaN compares as equal to itself.

.. function:: int fmpr_cmp(const fmpr_t x, const fmpr_t y)

    Returns negative, zero, or positive, depending on whether *x* is
    respectively smaller, equal, or greater compared to *y*.
    Comparison with NaN is undefined.

.. function:: int fmpr_cmpabs(const fmpr_t x, const fmpr_t y)

.. function:: int fmpr_cmpabs_ui(const fmpr_t x, ulong y)

    Compares the absolute values of *x* and *y*.

.. function:: int fmpr_cmp_2exp_si(const fmpr_t x, slong e)

.. function:: int fmpr_cmpabs_2exp_si(const fmpr_t x, slong e)

    Compares *x* (respectively its absolute value) with `2^e`.

.. function:: int fmpr_sgn(const fmpr_t x)

    Returns `-1`, `0` or `+1` according to the sign of *x*. The sign
    of NaN is undefined.

.. function:: void fmpr_min(fmpr_t z, const fmpr_t a, const fmpr_t b)

.. function:: void fmpr_max(fmpr_t z, const fmpr_t a, const fmpr_t b)

    Sets *z* respectively to the minimum and the maximum of *a* and *b*.

.. function:: slong fmpr_bits(const fmpr_t x)

    Returns the number of bits needed to represent the absolute value
    of the mantissa of *x*, i.e. the minimum precision sufficient to represent
    *x* exactly. Returns 0 if *x* is a special value.

.. function:: int fmpr_is_int(const fmpr_t x)

    Returns nonzero iff *x* is integer-valued.

.. function:: int fmpr_is_int_2exp_si(const fmpr_t x, slong e)

    Returns nonzero iff *x* equals `n 2^e` for some integer *n*.


Random number generation
-------------------------------------------------------------------------------

.. function:: void fmpr_randtest(fmpr_t x, flint_rand_t state, slong bits, slong mag_bits)

    Generates a finite random number whose mantissa has precision at most
    *bits* and whose exponent has at most *mag_bits* bits. The
    values are distributed non-uniformly: special bit patterns are generated
    with high probability in order to allow the test code to exercise corner
    cases.

.. function:: void fmpr_randtest_not_zero(fmpr_t x, flint_rand_t state, slong bits, slong mag_bits)

    Identical to *fmpr_randtest*, except that zero is never produced
    as an output.

.. function:: void fmpr_randtest_special(fmpr_t x, flint_rand_t state, slong bits, slong mag_bits)

    Identical to *fmpr_randtest*, except that the output occasionally
    is set to an infinity or NaN.


Input and output
-------------------------------------------------------------------------------

.. function:: void fmpr_print(const fmpr_t x)

    Prints the mantissa and exponent of *x* as integers, precisely showing
    the internal representation.

.. function:: void fmpr_printd(const fmpr_t x, slong digits)

    Prints *x* as a decimal floating-point number, rounding to the specified
    number of digits. This function is currently implemented using MPFR,
    and does not support large exponents.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void fmpr_neg(fmpr_t y, const fmpr_t x)

    Sets *y* to the negation of *x*.

.. function:: slong fmpr_neg_round(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)

    Sets *y* to the negation of *x*, rounding the result.

.. function:: void fmpr_abs(fmpr_t y, const fmpr_t x)

    Sets *y* to the absolute value of *x*.

.. function:: slong fmpr_add(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_add_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_add_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_add_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd)

    Sets `z = x + y`, rounded according to *prec* and *rnd*. The precision
    can be *FMPR_PREC_EXACT* to perform an exact addition, provided that the
    result fits in memory.

.. function:: slong _fmpr_add_eps(fmpr_t z, const fmpr_t x, int sign, slong prec, fmpr_rnd_t rnd)

    Sets *z* to the value that results by adding an infinitesimal quantity
    of the given sign to *x*, and rounding. The result is undefined
    if *x* is zero.

.. function:: slong fmpr_sub(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_sub_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_sub_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_sub_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd)

    Sets `z = x - y`, rounded according to *prec* and *rnd*. The precision
    can be  *FMPR_PREC_EXACT* to perform an exact addition, provided that the
    result fits in memory.

.. function:: slong fmpr_mul(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_mul_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_mul_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_mul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd)

    Sets `z = x \times y`, rounded according to prec and rnd. The precision
    can be *FMPR_PREC_EXACT* to perform an exact multiplication, provided that the
    result fits in memory.

.. function:: void fmpr_mul_2exp_si(fmpr_t y, const fmpr_t x, slong e)

.. function:: void fmpr_mul_2exp_fmpz(fmpr_t y, const fmpr_t x, const fmpz_t e)

    Sets *y* to *x* multiplied by `2^e` without rounding.

.. function:: slong fmpr_div(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_div_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_ui_div(fmpr_t z, ulong x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_div_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_si_div(fmpr_t z, slong x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_div_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_fmpz_div(fmpr_t z, const fmpz_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_fmpz_div_fmpz(fmpr_t z, const fmpz_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd)

    Sets `z = x / y`, rounded according to *prec* and *rnd*. If *y* is zero,
    *z* is set to NaN.

.. function:: slong fmpr_addmul(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_addmul_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_addmul_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_addmul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd)

    Sets `z = z + x \times y`, rounded according to *prec* and *rnd*. The
    intermediate multiplication is always performed without roundoff. The
    precision can be *FMPR_PREC_EXACT* to perform an exact addition, provided
    that the result fits in memory.

.. function:: slong fmpr_submul(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_submul_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_submul_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd)

.. function:: slong fmpr_submul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd)

    Sets `z = z - x \times y`, rounded according to *prec* and *rnd*. The
    intermediate multiplication is always performed without roundoff. The
    precision can be *FMPR_PREC_EXACT* to perform an exact subtraction, provided
    that the result fits in memory.

.. function:: slong fmpr_sqrt(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)

    Sets *z* to the square root of *x*, rounded according to *prec* and *rnd*.
    The result is NaN if *x* is negative.

.. function:: slong fmpr_rsqrt(fmpr_t z, const fmpr_t x, slong prec, fmpr_rnd_t rnd)

    Sets *z* to the reciprocal square root of *x*, rounded according to
    *prec* and *rnd*. The result is NaN if *x* is negative.
    At high precision, this is faster than computing a square root.

.. function:: slong fmpr_root(fmpr_t z, const fmpr_t x, ulong k, slong prec, fmpr_rnd_t rnd)

    Sets *z* to the *k*-th root of *x*, rounded to *prec* bits
    in the direction *rnd*.
    Warning: this function wraps MPFR, and is currently only fast for small *k*.

.. function:: void fmpr_pow_sloppy_fmpz(fmpr_t y, const fmpr_t b, const fmpz_t e, slong prec, fmpr_rnd_t rnd)

.. function:: void fmpr_pow_sloppy_ui(fmpr_t y, const fmpr_t b, ulong e, slong prec, fmpr_rnd_t rnd)

.. function:: void fmpr_pow_sloppy_si(fmpr_t y, const fmpr_t b, slong e, slong prec, fmpr_rnd_t rnd)

    Sets `y = b^e`, computed using without guaranteeing correct (optimal)
    rounding, but guaranteeing that the result is a correct upper or lower
    bound if the rounding is directional. Currently requires `b \ge 0`.


Special functions
-------------------------------------------------------------------------------


.. function:: slong fmpr_log(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)

    Sets *y* to `\log(x)`, rounded according to *prec* and *rnd*.
    The result is NaN if *x* is negative.
    This function is currently implemented using MPFR and does not
    support large exponents.

.. function:: slong fmpr_log1p(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)

    Sets *y* to `\log(1+x)`, rounded according to *prec* and *rnd*.
    This function
    computes an accurate value when *x* is small.
    The result is NaN if `1+x` is negative.
    This function is currently implemented using MPFR and does not
    support large exponents.

.. function:: slong fmpr_exp(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)

    Sets *y* to `\exp(x)`, rounded according to *prec* and *rnd*.
    This function is currently implemented using MPFR and does not
    support large exponents.

.. function:: slong fmpr_expm1(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)

    Sets *y* to `\exp(x)-1`, rounded according to *prec* and *rnd*.
    This function computes an accurate value when *x* is small.
    This function is currently implemented using MPFR and does not
    support large exponents.

