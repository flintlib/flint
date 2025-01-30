.. _mag:

**mag.h** -- fixed-precision unsigned floating-point numbers for bounds
===============================================================================

The :type:`mag_t` type holds an unsigned floating-point number with a
fixed-precision mantissa (30 bits) and an arbitrary-precision
exponent (represented as an :type:`fmpz_t`), suited for
representing magnitude bounds.
The special values zero and positive infinity are supported, but not NaN.

Operations that involve rounding will always produce a valid upper bound,
or a lower bound if the function name has the suffix *lower*.
For performance reasons, no attempt is made to compute the best possible bounds:
in general, a bound may be several ulps larger/smaller than the optimal bound.
Some functions such as :func:`mag_set` and :func:`mag_mul_2exp_si` are always
exact and therefore do not require separate *lower* versions.

A common mistake is to forget computing a lower bound for the argument
of a decreasing function that is meant to be bounded from above,
or vice versa. For example, to compute an upper bound for `(x+1)/(y+1)`,
the parameter *x* should initially be an upper bound while *y* should be
a lower bound, and one should do::

    mag_add_ui(tmp1, x, 1);
    mag_add_ui_lower(tmp2, y, 1);
    mag_div(res, tmp1, tmp2);

For a lower bound of the same expression, *x* should be a lower bound while
*y* should be an upper bound, and one should do::

    mag_add_ui_lower(tmp1, x, 1);
    mag_add_ui(tmp2, y, 1);
    mag_div_lower(res, tmp1, tmp2);

Applications requiring floating-point arithmetic with more flexibility
(such as correct rounding, or higher precision) should use the :type:`arf_t`
type instead. For calculations where a complex alternation between upper and
lower bounds is necessary, it may be cleaner to use :type:`arb_t`
arithmetic and convert to a :type:`mag_t` bound only in the end.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: mag_struct

    A :type:`mag_struct` holds a mantissa and an exponent.
    Special values are encoded by the mantissa being set to zero.

.. type:: mag_t

    A :type:`mag_t` is defined as an array of length one of type
    :type:`mag_struct`, permitting a :type:`mag_t` to be passed by reference.

Memory management
-------------------------------------------------------------------------------

.. function:: void mag_init(mag_t x)

    Initializes the variable *x* for use. Its value is set to zero.

.. function:: void mag_clear(mag_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

.. function:: void mag_swap(mag_t x, mag_t y)

    Swaps *x* and *y* efficiently.

.. function:: mag_ptr _mag_vec_init(slong n)

    Allocates a vector of length *n*. All entries are set to zero.

.. function:: void _mag_vec_clear(mag_ptr v, slong n)

    Clears a vector of length *n*.

.. function:: slong mag_allocated_bytes(const mag_t x)

    Returns the total number of bytes heap-allocated internally by this object.
    The count excludes the size of the structure itself. Add
    ``sizeof(mag_struct)`` to get the size of the object as a whole.

Special values
-------------------------------------------------------------------------------

.. function:: void mag_zero(mag_t res)

    Sets *res* to zero.

.. function:: void mag_one(mag_t res)

    Sets *res* to one.

.. function:: void mag_inf(mag_t res)

    Sets *res* to positive infinity.

.. function:: int mag_is_special(const mag_t x)

    Returns nonzero iff *x* is zero or positive infinity.

.. function:: int mag_is_zero(const mag_t x)

    Returns nonzero iff *x* is zero.

.. function:: int mag_is_inf(const mag_t x)

    Returns nonzero iff *x* is positive infinity.

.. function:: int mag_is_finite(const mag_t x)

    Returns nonzero iff *x* is not positive infinity (since there is no
    NaN value, this function is exactly the logical negation of :func:`mag_is_inf`).

Assignment and conversions
-------------------------------------------------------------------------------

.. function:: void mag_init_set(mag_t res, const mag_t x)

    Initializes *res* and sets it to the value of *x*. This operation is always exact.

.. function:: void mag_set(mag_t res, const mag_t x)

    Sets *res* to the value of *x*. This operation is always exact.

.. function:: void mag_set_d(mag_t res, double x)

.. function:: void mag_set_ui(mag_t res, ulong x)

.. function:: void mag_set_fmpz(mag_t res, const fmpz_t x)

    Sets *res* to an upper bound for `|x|`. The operation may be inexact
    even if *x* is exactly representable.

.. function:: void mag_set_d_lower(mag_t res, double x)

.. function:: void mag_set_ui_lower(mag_t res, ulong x)

.. function:: void mag_set_fmpz_lower(mag_t res, const fmpz_t x)

    Sets *res* to a lower bound for `|x|`.
    The operation may be inexact even if *x* is exactly representable.

.. function:: void mag_set_d_2exp_fmpz(mag_t res, double x, const fmpz_t y)

.. function:: void mag_set_fmpz_2exp_fmpz(mag_t res, const fmpz_t x, const fmpz_t y)

.. function:: void mag_set_ui_2exp_si(mag_t res, ulong x, slong y)

    Sets *res* to an upper bound for `|x| \cdot 2^y`.

.. function:: void mag_set_d_2exp_fmpz_lower(mag_t res, double x, const fmpz_t y)

.. function:: void mag_set_fmpz_2exp_fmpz_lower(mag_t res, const fmpz_t x, const fmpz_t y)

    Sets *res* to a lower bound for `|x| \cdot 2^y`.

.. function:: double mag_get_d(const mag_t x)

    Returns a *double* giving an upper bound for *x*.

.. function:: double mag_get_d_log2_approx(const mag_t x)

    Returns a *double* approximating `\log_2(x)`, suitable for estimating
    magnitudes (warning: not a rigorous bound).
    The value is clamped between *COEFF_MIN* and *COEFF_MAX*.

.. function:: void mag_get_fmpq(fmpq_t res, const mag_t x)

.. function:: void mag_get_fmpz(fmpz_t res, const mag_t x)

.. function:: void mag_get_fmpz_lower(fmpz_t res, const mag_t x)

    Sets *res*, respectively, to the exact rational number represented by *x*,
    the integer exactly representing the ceiling function of *x*, or the
    integer exactly representing the floor function of *x*.

    These functions are unsafe: the user must check in advance that *x* is of
    reasonable magnitude. If *x* is infinite or has a bignum exponent, an
    abort will be raised. If the exponent otherwise is too large or too small,
    the available memory could be exhausted resulting in undefined behavior.

Comparisons
-------------------------------------------------------------------------------

.. function:: int mag_equal(const mag_t x, const mag_t y)

    Returns nonzero iff *x* and *y* have the same value.

.. function:: int mag_cmp(const mag_t x, const mag_t y)

    Returns negative, zero, or positive, depending on whether *x*
    is smaller, equal, or larger than *y*.

.. function:: int mag_cmp_2exp_si(const mag_t x, slong y)

    Returns negative, zero, or positive, depending on whether *x*
    is smaller, equal, or larger than `2^y`.

.. function:: void mag_min(mag_t res, const mag_t x, const mag_t y)

.. function:: void mag_max(mag_t res, const mag_t x, const mag_t y)

    Sets *res* respectively to the smaller or the larger of *x* and *y*.

Input and output
-------------------------------------------------------------------------------

.. function:: void mag_print(const mag_t x)

    Prints *x* to standard output.

.. function:: void mag_fprint(FILE * file, const mag_t x)

    Prints *x* to the stream *file*.

.. function:: char * mag_dump_str(const mag_t x)

    Allocates a string and writes a binary representation of *x* to it that can
    be read by :func:`mag_load_str`. The returned string needs to be
    deallocated with *flint_free*.

.. function:: int mag_load_str(mag_t x, const char * str)

    Parses *str* into *x*. Returns a nonzero value if *str* is not formatted
    correctly.

.. function:: int mag_dump_file(FILE * stream, const mag_t x)

    Writes a binary representation of *x* to *stream* that can be read by
    :func:`mag_load_file`. Returns a nonzero value if the data could not be
    written.

.. function:: int mag_load_file(mag_t x, FILE * stream)

    Reads *x* from *stream*. Returns a nonzero value if the data is not
    formatted correctly or the read failed. Note that the data is assumed to be
    delimited by a whitespace or end-of-file, i.e., when writing multiple
    values with :func:`mag_dump_file` make sure to insert a whitespace to
    separate consecutive values.

Random generation
-------------------------------------------------------------------------------

.. function:: void mag_randtest(mag_t res, flint_rand_t state, slong expbits)

    Sets *res* to a random finite value, with an exponent up to *expbits* bits large.

.. function:: void mag_randtest_special(mag_t res, flint_rand_t state, slong expbits)

    Like :func:`mag_randtest`, but also sometimes sets *res* to infinity.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void mag_add(mag_t res, const mag_t x, const mag_t y)

.. function:: void mag_add_ui(mag_t res, const mag_t x, ulong y)

    Sets *res* to an upper bound for `x + y`.

.. function:: void mag_add_lower(mag_t res, const mag_t x, const mag_t y)

.. function:: void mag_add_ui_lower(mag_t res, const mag_t x, ulong y)

    Sets *res* to a lower bound for `x + y`.

.. function:: void mag_add_2exp_fmpz(mag_t res, const mag_t x, const fmpz_t e)

    Sets *res* to an upper bound for `x + 2^e`.

.. function:: void mag_add_ui_2exp_si(mag_t res, const mag_t x, ulong y, slong e)

    Sets *res* to an upper bound for `x + y 2^e`.

.. function:: void mag_sub(mag_t res, const mag_t x, const mag_t y)

    Sets *res* to an upper bound for `\max(x-y, 0)`.

.. function:: void mag_sub_lower(mag_t res, const mag_t x, const mag_t y)

    Sets *res* to a lower bound for `\max(x-y, 0)`.

.. function:: void mag_mul_2exp_si(mag_t res, const mag_t x, slong y)

.. function:: void mag_mul_2exp_fmpz(mag_t res, const mag_t x, const fmpz_t y)

    Sets *res* to `x \cdot 2^y`. This operation is exact.

.. function:: void mag_mul(mag_t res, const mag_t x, const mag_t y)

.. function:: void mag_mul_ui(mag_t res, const mag_t x, ulong y)

.. function:: void mag_mul_fmpz(mag_t res, const mag_t x, const fmpz_t y)

    Sets *res* to an upper bound for `xy`.

.. function:: void mag_mul_lower(mag_t res, const mag_t x, const mag_t y)

.. function:: void mag_mul_ui_lower(mag_t res, const mag_t x, ulong y)

.. function:: void mag_mul_fmpz_lower(mag_t res, const mag_t x, const fmpz_t y)

    Sets *res* to a lower bound for `xy`.

.. function:: void mag_addmul(mag_t z, const mag_t x, const mag_t y)

    Sets *z* to an upper bound for `z + xy`.

.. function:: void mag_div(mag_t res, const mag_t x, const mag_t y)

.. function:: void mag_div_ui(mag_t res, const mag_t x, ulong y)

.. function:: void mag_div_fmpz(mag_t res, const mag_t x, const fmpz_t y)

    Sets *res* to an upper bound for `x / y`.

.. function:: void mag_div_lower(mag_t res, const mag_t x, const mag_t y)

    Sets *res* to a lower bound for `x / y`.

.. function:: void mag_inv(mag_t res, const mag_t x)

    Sets *res* to an upper bound for `1 / x`.

.. function:: void mag_inv_lower(mag_t res, const mag_t x)

    Sets *res* to a lower bound for `1 / x`.


Fast, unsafe arithmetic
-------------------------------------------------------------------------------

The following methods assume that all inputs are finite and that all exponents
(in all inputs as well as the final result) fit as *fmpz* inline values.
They also assume that the output variables do not have promoted exponents,
as they will be overwritten directly (thus leaking memory).

.. function:: void mag_fast_init_set(mag_t x, const mag_t y)

    Initialises *x* and sets it to the value of *y*.

.. function:: void mag_fast_zero(mag_t res)

    Sets *res* to zero.

.. function:: int mag_fast_is_zero(const mag_t x)

    Returns nonzero iff *x* to zero.

.. function:: void mag_fast_mul(mag_t res, const mag_t x, const mag_t y)

    Sets *res* to an upper bound for `xy`.

.. function:: void mag_fast_addmul(mag_t z, const mag_t x, const mag_t y)

    Sets *z* to an upper bound for `z + xy`.

.. function:: void mag_fast_add_2exp_si(mag_t res, const mag_t x, slong e)

    Sets *res* to an upper bound for `x + 2^e`.

.. function:: void mag_fast_mul_2exp_si(mag_t res, const mag_t x, slong e)

    Sets *res* to an upper bound for `x 2^e`.

Powers and logarithms
-------------------------------------------------------------------------------

.. function:: void mag_pow_ui(mag_t res, const mag_t x, ulong e)

.. function:: void mag_pow_fmpz(mag_t res, const mag_t x, const fmpz_t e)

    Sets *res* to an upper bound for `x^e`.

.. function:: void mag_pow_ui_lower(mag_t res, const mag_t x, ulong e)

.. function:: void mag_pow_fmpz_lower(mag_t res, const mag_t x, const fmpz_t e)

    Sets *res* to a lower bound for `x^e`.

.. function:: void mag_sqrt(mag_t res, const mag_t x)

    Sets *res* to an upper bound for `\sqrt{x}`.

.. function:: void mag_sqrt_lower(mag_t res, const mag_t x)

    Sets *res* to a lower bound for `\sqrt{x}`.

.. function:: void mag_rsqrt(mag_t res, const mag_t x)

    Sets *res* to an upper bound for `1/\sqrt{x}`.

.. function:: void mag_rsqrt_lower(mag_t res, const mag_t x)

    Sets *res* to an lower bound for `1/\sqrt{x}`.

.. function:: void mag_hypot(mag_t res, const mag_t x, const mag_t y)

    Sets *res* to an upper bound for `\sqrt{x^2 + y^2}`.

.. function:: void mag_root(mag_t res, const mag_t x, ulong n)

    Sets *res* to an upper bound for `x^{1/n}`. 

.. function:: void mag_log(mag_t res, const mag_t x)

    Sets *res* to an upper bound for `\log(\max(1,x))`.

.. function:: void mag_log_lower(mag_t res, const mag_t x)

    Sets *res* to a lower bound for `\log(\max(1,x))`.

.. function:: void mag_neg_log(mag_t res, const mag_t x)

    Sets *res* to an upper bound for `-\log(\min(1,x))`, i.e. an upper
    bound for `|\log(x)|` for `x \le 1`.

.. function:: void mag_neg_log_lower(mag_t res, const mag_t x)

    Sets *res* to a lower bound for `-\log(\min(1,x))`, i.e. a lower
    bound for `|\log(x)|` for `x \le 1`.

.. function:: void mag_log_ui(mag_t res, ulong n)

    Sets *res* to an upper bound for `\log(n)`.

.. function:: void mag_log1p(mag_t res, const mag_t x)

    Sets *res* to an upper bound for `\log(1+x)`. The bound is computed
    accurately for small *x*.

.. function:: void mag_exp(mag_t res, const mag_t x)

    Sets *res* to an upper bound for `\exp(x)`.

.. function:: void mag_exp_lower(mag_t res, const mag_t x)

    Sets *res* to a lower bound for `\exp(x)`.

.. function:: void mag_expinv(mag_t res, const mag_t x)

    Sets *res* to an upper bound for `\exp(-x)`.

.. function:: void mag_expinv_lower(mag_t res, const mag_t x)

    Sets *res* to a lower bound for `\exp(-x)`.

.. function:: void mag_expm1(mag_t res, const mag_t x)

    Sets *res* to an upper bound for `\exp(x) - 1`. The bound is computed
    accurately for small *x*.

.. function:: void mag_exp_tail(mag_t res, const mag_t x, ulong N)

    Sets *res* to an upper bound for `\sum_{k=N}^{\infty} x^k / k!`.

.. function:: void mag_binpow_uiui(mag_t res, ulong m, ulong n)

    Sets *res* to an upper bound for `(1 + 1/m)^n`.

.. function:: void mag_geom_series(mag_t res, const mag_t x, ulong N)

    Sets *res* to an upper bound for `\sum_{k=N}^{\infty} x^k`.

Special functions
-------------------------------------------------------------------------------

.. function:: void mag_const_pi(mag_t res)

.. function:: void mag_const_pi_lower(mag_t res)

    Sets *res* to an upper (respectively lower) bound for `\pi`.

.. function:: void mag_atan(mag_t res, const mag_t x)

.. function:: void mag_atan_lower(mag_t res, const mag_t x)

    Sets *res* to an upper (respectively lower) bound for `\operatorname{atan}(x)`.

.. function:: void mag_cosh(mag_t res, const mag_t x)

.. function:: void mag_cosh_lower(mag_t res, const mag_t x)

.. function:: void mag_sinh(mag_t res, const mag_t x)

.. function:: void mag_sinh_lower(mag_t res, const mag_t x)

    Sets *res* to an upper or lower bound for `\cosh(x)` or `\sinh(x)`.

.. function:: void mag_fac_ui(mag_t res, ulong n)

    Sets *res* to an upper bound for `n!`.

.. function:: void mag_rfac_ui(mag_t res, ulong n)

    Sets *res* to an upper bound for `1/n!`.

.. function:: void mag_bin_uiui(mag_t res, ulong n, ulong k)

    Sets *res* to an upper bound for the binomial coefficient `{n \choose k}`.

.. function:: void mag_bernoulli_div_fac_ui(mag_t res, ulong n)

    Sets *res* to an upper bound for `|B_n| / n!` where `B_n` denotes
    a Bernoulli number.

.. function:: void mag_polylog_tail(mag_t res, const mag_t z, slong s, ulong d, ulong N)

    Sets *res* to an upper bound for

    .. math::

        \sum_{k=N}^{\infty} \frac{z^k \log^d(k)}{k^s}.

    The bounding strategy is described in :ref:`algorithms_polylogarithms`.
    Note: in applications where `s` in this formula may be
    real or complex, the user can simply
    substitute any convenient integer `s'` such that `s' \le \operatorname{Re}(s)`.

.. function:: void mag_hurwitz_zeta_uiui(mag_t res, ulong s, ulong a)

    Sets *res* to an upper bound for `\zeta(s,a) = \sum_{k=0}^{\infty} (k+a)^{-s}`.
    We use the formula

    .. math::

        \zeta(s,a) \le \frac{1}{a^s} + \frac{1}{(s-1) a^{s-1}}

    which is obtained by estimating the sum by an integral.
    If `s \le 1` or `a = 0`, the bound is infinite.

