.. _nfloat:

**nfloat.h** -- packed floating-point numbers with n-word precision
===============================================================================

This module provides binary floating-point numbers in a flat representation
with precision in small fixed multiples of the word size
(64, 128, 192, 256, ... bits on a 64-bit machine). The exponent range is
close to a full word.
A number with `n`-limb precision is stored as `n+2` contiguous
limbs as follows:

    +---------------+
    | exponent limb |
    +---------------+
    |   sign limb   |
    +---------------+
    |  mantissa[0]  |
    +---------------+
    |      ...      |
    +---------------+
    | mantissa[n-1] |
    +---------------+

For normal (nonzero and finite) values `x`,
the most significant limb of the mantissa is always normalised
to have its most significant bit set, and the exponent `e` is the
unique integer such that `|x| \in [0.5, 1) \cdot 2^e`.
Special (zero or nonfinite) values are encoded using special
values of the exponent field, with junk data in the mantissa.

This type has the advantage that floating-point numbers with the
same precision can be packed together tightly in vectors
and created on the stack without heap allocation.
The precision of an ``nfloat`` context object and its elements cannot
be changed; to switch precision, one must convert to a different
context object. For higher precision than supported by ``nfloat``
and for calculations that require fine-grained precision adjustments,
one should use :type:`arf_t` instead.

The focus is on fast calculation, not bitwise-defined results.
Without directed rounding, atomic operations give the correctly rounded
(truncated) result or slightly worse, e.g. with 1-2 ulp error. With the
directed rounding flags, most atomic operations are correctly rounded,
and the other supported operations return valid bounds; see
:ref:`nfloat-rounding` for the details of each operation. Results
that are not correctly rounded are not guaranteed to be identical on
64-bit and 32-bit machines.
Planned features include:

* Support for special values (partially implemented)
* Optional (but slower) IEEE 754-like semantics
* Complex and ball types

This module is designed to use the :ref:`generics <gr>` interface.
As such, the domain is represented by a :type:`gr_ctx_t` context object,
methods return status flags (``GR_SUCCESS``, ``GR_UNABLE``, ``GR_DOMAIN``),
and one can use generic structures such as :type:`gr_poly_t` for
polynomials and :type:`gr_mat_t` for matrices.

Types, macros and constants
-------------------------------------------------------------------------------

.. macro :: NFLOAT_MIN_LIMBS
            NFLOAT_MAX_LIMBS

    The number of limbs `n` permitted as precision. The current
    limits are are `1 \le n \le 66` on a 64-bit machine and
    `1 \le n \le 132` on a 32-bit machine, permitting precision
    up to 4224 bits. The upper limit exists so that elements and
    temporary buffers are safe to allocate on the stack and so that
    simple operations like swapping are not too expensive.

.. type:: nfloat_ptr
          nfloat_srcptr

    Pointer to an ``nfloat`` element or vector of elements of any
    precision. Since this is a void type, one must cast to the correct
    size before doing pointer arithmetic, e.g. via the ``GR_ENTRY``
    macro.

.. type:: nfloat64_struct
          nfloat128_struct
          nfloat192_struct
          nfloat256_struct
          nfloat384_struct
          nfloat512_struct
          nfloat1024_struct
          nfloat2048_struct
          nfloat4096_struct
          nfloat64_t
          nfloat128_t
          nfloat192_t
          nfloat256_t
          nfloat384_t
          nfloat512_t
          nfloat1024_t
          nfloat2048_t
          nfloat4096_t

    For convenience we define types of the correct structure size for
    some common levels of bit precision. An ``nfloatX_t`` is defined as
    a length-one array of ``nfloatX_struct``, permitting it to be
    passed by reference.

    Sample usage:

    .. code-block:: c

        gr_ctx_t ctx;
        nfloat256_t x, y;

        nfloat_ctx_init(ctx, 256, 0);   /* precision must match the type */
        gr_init(x, ctx);
        gr_init(y, ctx);

        gr_ctx_println(ctx);

        GR_MUST_SUCCEED(gr_set_ui(x, 5, ctx));
        GR_MUST_SUCCEED(gr_set_ui(y, 7, ctx));
        GR_MUST_SUCCEED(gr_div(x, x, y, ctx));
        GR_MUST_SUCCEED(gr_println(x, ctx));

        gr_clear(x, ctx);
        gr_clear(y, ctx);
        gr_ctx_clear(ctx);

.. macro:: NFLOAT_HEADER_LIMBS
           NFLOAT_EXP(x)
           NFLOAT_SGNBIT(x)
           NFLOAT_D(x)
           NFLOAT_DATA(x)

.. macro:: NFLOAT_MAX_ALLOC

.. macro:: NFLOAT_MIN_EXP
           NFLOAT_MAX_EXP

.. macro:: NFLOAT_EXP_ZERO
           NFLOAT_EXP_POS_INF
           NFLOAT_EXP_NEG_INF
           NFLOAT_EXP_NAN
           NFLOAT_IS_SPECIAL(x)
           NFLOAT_IS_ZERO(x)
           NFLOAT_IS_POS_INF(x)
           NFLOAT_IS_NEG_INF(x)
           NFLOAT_IS_INF(x)
           NFLOAT_IS_NAN(x)

Context objects
-------------------------------------------------------------------------------

.. function:: int nfloat_ctx_init(gr_ctx_t ctx, slong prec, int flags)

    Initializes *ctx* to represent a domain of floating-point numbers
    with bit precision *prec* rounded up to a full word
    (for example, ``prec = 53`` actually creates a domain with
    64-bit precision).

    Returns ``GR_UNABLE`` without initializing the context object
    if the given precision is too large to be supported, otherwise
    returns ``GR_SUCCESS``.

    Admissible flags are listed below.

.. macro:: NFLOAT_ALLOW_UNDERFLOW

    By default, operations that would underflow the exponent range
    output a garbage value and return ``GR_UNABLE``.
    Setting this flag allows such operations to
    output zero and return ``GR_SUCCESS`` instead.

.. macro:: NFLOAT_ALLOW_INF

    Allow creation of infinities.
    By default, operations that would overflow the exponent range
    output a garbage value and return ``GR_UNABLE`` or ``GR_DOMAIN``.
    Setting this flag allows such operations to
    output an infinity and return ``GR_SUCCESS`` instead.

.. macro:: NFLOAT_ALLOW_NAN

    Allow creation of NaNs.
    By default, operations that are meaningless
    output a garbage value and return ``GR_UNABLE`` or ``GR_DOMAIN``.
    Setting this flag allows such operations to
    output NaN and return ``GR_SUCCESS`` instead.

.. macro:: NFLOAT_RND_FLOOR
           NFLOAT_RND_CEIL

    Enable directed rounding towards `-\infty` or `+\infty` respectively.
    At most one of these flags may be set. If neither is set, rounding is done
    in an arbitrary direction for efficiency
    (usually but not always truncating towards zero).
    See below for a list of operations which support directed
    rounding.

Infinities and NaNs are disabled by default to improve performance,
as this allows certain functions to skip checks for such values.

.. function:: int nfloat_ctx_set_func_prec(gr_ctx_t ctx, slong prec)
              slong nfloat_ctx_get_func_prec(gr_ctx_t ctx)

    Sets (gets) the *function precision* of *ctx*: the precision in bits
    to which the elementary functions (exponentials, logarithms,
    trigonometric and hyperbolic functions and their inverses, powers)
    are computed. It defaults to the full precision ``FLINT_BITS *
    nlimbs`` of the context; a smaller value is clamped to at least 1 bit
    and makes these functions cheaper. For example, with ``nfloat64``
    and a function precision of 50 bits, most functions evaluate their
    kernels with one limb instead of two. The precision of the
    arithmetic operations is unaffected. Setting a function precision
    less than 1 returns ``GR_UNABLE``.

Basic operations and arithmetic
-------------------------------------------------------------------------------

Basic functionality for the ``gr`` method table.
These methods are interchangeable with their ``gr`` counterparts.

.. function:: int nfloat_ctx_write(gr_stream_t out, gr_ctx_t ctx)

.. function:: void nfloat_init(nfloat_ptr res, gr_ctx_t ctx)

    Initializes *res* to the zero element.

.. function:: void nfloat_clear(nfloat_ptr res, gr_ctx_t ctx)

    Since ``nfloat`` elements do no allocation, this is a no-op.

.. function:: void nfloat_swap(nfloat_ptr x, nfloat_ptr y, gr_ctx_t ctx)

.. function:: int nfloat_set(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

.. function:: truth_t nfloat_equal(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

.. function:: int nfloat_ctx_set_real_prec(gr_ctx_t ctx, slong prec)

    Since ``nfloat`` contexts do not allow variable precision,
    this does nothing and returns ``GR_UNABLE``.

.. function:: int nfloat_ctx_get_real_prec(slong * res, gr_ctx_t ctx)

    Sets *res* to the precision in bits and returns ``GR_SUCCESS``.

.. function:: int nfloat_zero(nfloat_ptr res, gr_ctx_t ctx)
              int nfloat_one(nfloat_ptr res, gr_ctx_t ctx)
              int nfloat_neg_one(nfloat_ptr res, gr_ctx_t ctx)
              int nfloat_pos_inf(nfloat_ptr res, gr_ctx_t ctx)
              int nfloat_neg_inf(nfloat_ptr res, gr_ctx_t ctx)
              int nfloat_nan(nfloat_ptr res, gr_ctx_t ctx)

.. function:: truth_t nfloat_is_zero(nfloat_srcptr x, gr_ctx_t ctx)
              truth_t nfloat_is_one(nfloat_srcptr x, gr_ctx_t ctx)
              truth_t nfloat_is_neg_one(nfloat_srcptr x, gr_ctx_t ctx)

.. function:: int nfloat_set_ui(nfloat_ptr res, ulong x, gr_ctx_t ctx)
              int nfloat_set_si(nfloat_ptr res, slong x, gr_ctx_t ctx)
              int nfloat_set_fmpz(nfloat_ptr res, const fmpz_t x, gr_ctx_t ctx)

.. function:: int _nfloat_set_mpn_2exp(nfloat_ptr res, nn_srcptr x, slong xn, slong exp, int xsgnbit, gr_ctx_t ctx)
              int nfloat_set_mpn_2exp(nfloat_ptr res, nn_srcptr x, slong xn, slong exp, int xsgnbit, gr_ctx_t ctx)

.. function:: int nfloat_set_arf(nfloat_ptr res, const arf_t x, gr_ctx_t ctx)
              int nfloat_get_arf(arf_t res, nfloat_srcptr x, gr_ctx_t ctx)

.. function:: int nfloat_get_d_2exp_si(double * m, slong * e, nfloat_srcptr x, gr_ctx_t ctx)

.. function:: int nfloat_get_fmpz(fmpz_t res, nfloat_srcptr x, gr_ctx_t ctx)

    Succeeds only if ``x`` is integer-valued.

.. function:: int nfloat_set_fmpq(nfloat_ptr res, const fmpq_t v, gr_ctx_t ctx)
              int nfloat_set_d(nfloat_ptr res, double x, gr_ctx_t ctx)
              int nfloat_set_str(nfloat_ptr res, const char * x, gr_ctx_t ctx)
              int nfloat_set_other(nfloat_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)

.. function:: int nfloat_write(gr_stream_t out, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_randtest(nfloat_ptr res, flint_rand_t state, gr_ctx_t ctx)

.. function:: int nfloat_cmp(int * res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
              int nfloat_cmpabs(int * res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

.. function:: int nfloat_neg(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_abs(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_add(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
              int nfloat_sub(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
              int nfloat_mul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
              int nfloat_submul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
              int nfloat_addmul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
              int nfloat_sqr(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

.. function:: int nfloat_mul_2exp_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx)

.. function:: int nfloat_inv(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
              int nfloat_div_ui(nfloat_ptr res, nfloat_srcptr x, ulong y, gr_ctx_t ctx)
              int nfloat_div_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx)

    Division. Without directed rounding, the result is the quotient
    truncated towards zero or (for divisors of three or more nonzero
    limbs) possibly one ulp larger in magnitude, so that the error is
    less than 1 ulp. With directed rounding, the result is the correctly
    rounded floor or ceiling of the exact quotient.
    Division by zero gives NaN (``GR_UNABLE`` unless NaNs are allowed).

    For one and two limbs, the mantissas are divided by a fully inlined
    2/1 respectively two 3/2 divisions, with the dividend shifted right
    by one bit beforehand when it is larger than the divisor so that
    the quotient comes out normalized; the remainder gives the exact
    rounding information. These use hardware division where it is fast
    (``FLINT_PREINVERT_LIMB_USE_NATIVE``) and a precomputed inverse
    otherwise. Divisors with a single nonzero limb (for instance,
    integers) use a chain of 2/1 divisions, `O(n)`.
    Otherwise, a schoolbook approximate division
    (:func:`_flint_mpn_divapprox_basecase_preinv1`) computes the
    quotient or the quotient plus one at the target precision; with
    directed rounding, one guard limb is computed and the truncation is
    certified unless the guard limb is 0 or 1 (in which case, for
    instance for exact quotients, the exact quotient and remainder are
    computed).

.. function:: int nfloat_sqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_rsqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

    Square root and reciprocal square root. With directed rounding, the
    result is the correctly rounded floor or ceiling of the exact value.

    The square root is computed exactly (truncated, with an exactness
    flag) by :func:`flint_mpn_sqrtrem` applied to the mantissa padded
    to `2n` limbs, so without directed rounding the result is the square
    root truncated towards zero. At one limb on 64-bit machines, the
    integer square root of the padded mantissa is instead computed in
    registers: a double-precision root (within `2^{12}` of the exact
    root), one correction from the exact residual evaluated in double
    precision (which leaves the floor or one off), and a final exact
    residual, which gives the floor and the exactness flag. At two limbs
    on 64-bit machines, the root is first approximated at three limbs
    from the two-limb reciprocal square root described below and one
    Newton step with the exact residual (within 3 units of the guard
    limb), which certifies the truncation (and inexactness) unless the
    guard limb is within that distance of 0 or `2^{64}`.

    Without directed rounding, the reciprocal square root has an error
    less than `1 + 2^{-37}` ulp. For one limb, a double-precision
    approximation is refined by one Newton step at 128-bit precision, with
    the residual evaluated exactly. At two limbs on 64-bit machines, the
    two-limb starting value below (within 8 units of its last limb) is
    refined in registers by one Newton step whose residual is evaluated
    exactly from the six-limb product of the mantissa and the square of
    the starting value, giving three limbs within 2 units of the last
    one. For `n \ge 3` limbs, a fixed-point Newton iteration (from a
    two-limb starting value computed in the same way, with an added
    third-order term) yields the result with one guard limb and an error
    bound of about one unit in the guard limb. With
    directed rounding, the truncation is certified using the error bounds,
    and in the rare uncertain cases the exact result
    `\lfloor \sqrt{\lfloor z \rfloor} \rfloor`, `z = 2^{3nw}/(4a)` or
    `2^{3nw}/(2a)` in terms of the `n`-limb mantissa `a` and the word size
    `w`, is computed.

.. function:: int nfloat_sgn(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_im(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

.. function:: int nfloat_floor(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_ceil(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_trunc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_nint(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

.. function:: int nfloat_pi(nfloat_ptr res, gr_ctx_t ctx)

Elementary functions
-------------------------------------------------------------------------------

The elementary functions are computed natively, using the fixed-point
kernels of the :type:`mp_real_t` module on `[0, 1)` (the bitwise and
reduced-argument kernels; ``nfloat`` precisions never call for bit-burst
evaluation). The argument reductions and the handling of tiny and huge
arguments work directly on the ``nfloat`` mantissas, in registers for
small sizes.

*Accuracy.* The functions are not correctly rounded. Without directed
rounding, the error is at most about one ulp at the function precision
of the context (see :func:`nfloat_ctx_set_func_prec`); with
``NFLOAT_RND_FLOOR`` or ``NFLOAT_RND_CEIL``, the result is a valid lower
or upper bound of the exact value (within a few ulp at the function
precision). Exact results are returned in exact special cases (for
example ``exp(0)``, ``log(1)``, ``sin_pi`` of multiples of 1/2, ``log2``
of powers of two). Results that over- or underflow the exponent range
follow the flags of the context. Arguments for which the functions are
undefined in the reals (for example ``log`` of a negative number, ``asin``
of `|x| > 1`) give NaN (``GR_UNABLE`` or ``GR_DOMAIN`` unless NaNs are
allowed).

On 64-bit machines, the functions :func:`nfloat_exp`, :func:`nfloat_expm1`,
:func:`nfloat_log`, :func:`nfloat_log1p`, :func:`nfloat_pow`,
:func:`nfloat_sin_cos`, :func:`nfloat_sin`, :func:`nfloat_cos`,
:func:`nfloat_tan`, :func:`nfloat_cot`, their versions at `\pi x`,
:func:`nfloat_sinh_cosh`, :func:`nfloat_sinh`, :func:`nfloat_cosh`,
:func:`nfloat_tanh`, :func:`nfloat_coth`, :func:`nfloat_atan` and
:func:`nfloat_atan2` at one and two limbs have in-register fast paths for
common arguments (for example every argument below the overflow
threshold for the exponential, `2^{-33} \le |x| < 2^{63}` for the sine
and cosine at one limb): the argument reduction on the mantissa, the
table-driven kernels of :type:`mp_real_t` at two resp. three limbs
(shared with its 128- and 192-bit functions) and the final
normalization, without temporary arrays. The hyperbolic functions
use a paired kernel for `\exp(t)` and `\exp(-t)` (and one division
for :func:`nfloat_tanh` and :func:`nfloat_coth`), the tangent and the
cotangent divide the sine and the cosine in registers, :func:`nfloat_pow` multiplies a logarithm with
one more limb (or two, for large `|y|`) by `y` exactly before the
exponential, and the arctangent of `|x| > 1` and :func:`nfloat_atan2`
evaluate the arctangent of a ratio of mantissas with a single division.
Arguments where the reduced argument would lose too much relative
accuracy (close to 1 for the logarithm, close to a multiple of `\pi/2`
for the sine, cosine and tangent, very small arguments of
:func:`nfloat_expm1`, :func:`nfloat_log1p` and :func:`nfloat_sinh`)
and the rarer ranges (for example :func:`nfloat_pow` with `|y| \ge 2^{40}`
at two limbs) take the general path.

From three to nine limbs the general path runs on the dense-table
kernels of :type:`mp_real_t` at four to ten limbs, and from ten to
thirty-one limbs on the dense-table kernels with tables computed at run time
(shared by all threads, built on first use; see the :type:`mp_real_t`
documentation); the arctangent of `|x| > 1` forms `1/|x|` by their
small reciprocal up to ten limbs.

Outside the fast paths, the functions wrap the ball functions of
:type:`mp_real_t` (for example :func:`mp_real_sin_cos_bits`,
:func:`mp_real_asinh_bits`), evaluated on the exact argument at
increasing precision until the ball determines the result: this
handles trigonometric functions of huge arguments (exponent above
65536, up to `2^{22}`, beyond which ``GR_UNABLE`` is returned), the
less common functions for arguments whose intermediate results over- or
underflow, and the maximum precision ``NFLOAT_MAX_LIMBS`` (which leaves
no room for a guard limb).

.. function:: int nfloat_exp(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_expm1(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_exp2(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_log(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_log1p(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_log2(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

    Exponentials and logarithms. The exponentials reduce `|x|` modulo
    `\log 2` in fixed point (with one integral limb); ``expm1`` and
    ``log1p`` have relative accuracy near zero, and ``exp2`` and ``log2``
    are exact for integer arguments respectively powers of two.

.. function:: int nfloat_sin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_cos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sin_cos(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_tan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sin_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_cos_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sin_cos_pi(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_tan_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

    Trigonometric functions. The argument is reduced modulo `\pi/2` in
    fixed point (repeated with more bits when the reduced argument is
    close to a zero of the output, so that the results have relative
    accuracy); the functions of `\pi x` reduce `x` modulo 2 exactly.

.. function:: int nfloat_atan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_atan2(nfloat_ptr res, nfloat_srcptr y, nfloat_srcptr x, gr_ctx_t ctx)

    Arctangents. For `|x| > 1`, ``atan`` uses `\pi/2 - \operatorname{atan}(1/|x|)`;
    ``atan2`` computes the ratio of the smaller and larger magnitude
    and adds the octant. ``atan2(0, 0)`` is defined as 0.

.. function:: int nfloat_sinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_cosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sinh_cosh(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_tanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

    Hyperbolic functions, from `E = \exp(|x|)` and `1/E` (with the
    difference `E - 1/E` in fixed point for `|x| < 1`).

.. function:: int nfloat_pow(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

    Power `x^y = \exp(y \log |x|)`, with the sign `(-1)^y` for `x < 0`
    and integer `y` (NaN for `x < 0` and noninteger `y`).
    Simple cases (`y = \pm 1, 2, \pm 1/2`, `|x| = 1`, integer powers of
    powers of two) are handled directly.

.. function:: int nfloat_exp10(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_log10(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_cot(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sec(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_csc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sinc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_cot_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sec_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_csc_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sinc_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_coth(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sech(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_csch(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_asin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_acos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_asin_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_acos_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_atan_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_acot(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_asec(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_acsc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_acot_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_asec_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_acsc_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_asinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_acosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_atanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_acoth(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_asech(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_acsch(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

    The less common functions, with compact code: after the exact
    special values (zeros, poles, `|x| = 1` for the inverse functions,
    multiples of 1/2 for the functions of `\pi x`), they are evaluated as
    short compositions of the functions above and arithmetic operations
    in an internal context with a guard limb (none when the function
    precision is sufficiently below the precision), using formulas
    without cancellation (for example `\operatorname{asin}(x) =
    \operatorname{atan2}(x, \sqrt{(1-x)(1+x)})` and
    `\operatorname{asinh}(x) = \operatorname{log1p}(|x| + x^2 / (1 + \sqrt{1 + x^2}))`)
    and a rigorous bound for the composed error, which is added to obtain
    valid bounds with directed rounding.

    On 64-bit machines, at one and two limbs, the reciprocal functions
    (``sec``, ``csc``, ``sech``, ``csch``, also of `\pi x`) and ``sinc``
    (also of `\pi x`) instead divide the fixed-point outputs of the
    in-register fast paths of the sine, cosine, and the paired
    exponentials, and the inverse functions are evaluated in fixed point
    at one more limb: `y^2 \pm 1` exactly, square roots by a refined
    reciprocal square root (as for :func:`nfloat_sqrt` at two limbs), the
    arctangent of a ratio (``asin``, ``acos``, ``asec``, ``acsc``,
    ``acot`` and the functions divided by `\pi`, which multiply by
    `2/\pi`) or the logarithm kernels (the inverse hyperbolic functions,
    as sums of logarithms of the same sign, for example
    `\operatorname{atanh}(y) = (\operatorname{log1p}(y) - \log(1 - y))/2`).
    Arguments where an intermediate would lose relative accuracy at that
    precision (results below about `2^{-32}` resp. `2^{-40}` from small
    differences, `|x| \ge 2^{62}`) take the general path.

.. function:: int nfloat_hypot(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

    Sets *res* to `\sqrt{x^2 + y^2}`, correctly rounded in every rounding
    mode (truncated towards zero without directed rounding), without
    intermediate overflow or underflow. With `|x| \ge |y|` and
    `d = e_x - e_y`, the result is `|x| (1 + \delta)` with
    `0 < \delta < 2^{-p-3}` (`p` the precision), rounded accordingly, when
    `d` exceeds about half the precision; otherwise `x^2 + y^2 = S 2^{2e}` for an integer `S` of at
    most `3n + 2` limbs (the squares of the mantissas, aligned exactly),
    whose truncated square root with an inexactness flag gives the result.
    At one limb on 64-bit machines, `S` is formed in registers and the
    truncated root is computed exactly from the integer square root of its
    top two limbs (as for :func:`nfloat_sqrt`), since the discarded low
    bits of `S` do not change the floor of the square root. At two limbs,
    `S` is formed in registers (the smaller square truncated to the four
    limbs of the larger) and its root is approximated with a guard limb,
    which certifies the truncation in all but rare cases. From 12
    limbs (or from 5 limbs when the exponents differ by at least `3/8` of
    the precision), the exact squares are replaced by high squarings
    (:func:`flint_mpn_sqrhigh`), whose sum is a lower bound within about
    `4n` units of its last limb, and the square root with a guard limb
    certifies the truncation likewise, the exact computation being the
    fallback. This also implements :func:`nfloat_complex_abs`.

.. function:: int nfloat_gamma(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_zeta(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

    These are computed via arb (with an arf result).

Vector functions
-------------------------------------------------------------------------------

Overrides for generic ``gr`` vector operations with inlined or partially inlined
code for reduced overhead.

.. function:: void _nfloat_vec_init(nfloat_ptr res, slong len, gr_ctx_t ctx)
              void _nfloat_vec_clear(nfloat_ptr res, slong len, gr_ctx_t ctx)
              int _nfloat_vec_set(nfloat_ptr res, nfloat_srcptr x, slong len, gr_ctx_t ctx)
              int _nfloat_vec_zero(nfloat_ptr res, slong len, gr_ctx_t ctx)

.. function:: int _nfloat_vec_add(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
              int _nfloat_vec_sub(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
              int _nfloat_vec_mul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
              int _nfloat_vec_mul_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)
              int _nfloat_vec_addmul_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)
              int _nfloat_vec_submul_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)

.. function:: int _nfloat_vec_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
              int _nfloat_vec_div_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr c, gr_ctx_t ctx)
              int _nfloat_vec_div_scalar_ui(nfloat_ptr res, nfloat_srcptr x, slong len, ulong c, gr_ctx_t ctx)
              int _nfloat_vec_div_scalar_si(nfloat_ptr res, nfloat_srcptr x, slong len, slong c, gr_ctx_t ctx)

    Vector division, with the same results as the corresponding scalar
    functions except that :func:`_nfloat_vec_div_scalar` for three or more
    limbs (and a divisor with two or more nonzero limbs) may give a
    different result within 1 ulp without directed rounding. The
    elementwise division inlines the one- and two-limb kernels. Division by
    a scalar precomputes the one-limb inverse, the 3/2 inverse, or for
    three or more limbs an approximate reciprocal with one guard limb,
    so that each quotient costs one high product of `n + 1` limbs; with
    directed rounding the truncation is certified using the guard limb and
    otherwise recomputed by exact division.

.. function:: int _nfloat_vec_dot(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
              int _nfloat_vec_dot_rev(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)

Matrix functions
-------------------------------------------------------------------------------

.. function:: int nfloat_mat_mul_fixed(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong max_extra_prec, gr_ctx_t ctx)
              int nfloat_mat_mul_block(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong min_block_size, gr_ctx_t ctx)
              int nfloat_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

    Different implementations of matrix multiplication.

.. function:: int nfloat_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int nfloat_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int nfloat_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
              int nfloat_mat_lq(gr_mat_t L, gr_mat_t Q, const gr_mat_t A, gr_ctx_t ctx)

Internal functions
-------------------------------------------------------------------------------

.. function:: int _nfloat_underflow(nfloat_ptr res, int sgnbit, gr_ctx_t ctx)
              int _nfloat_overflow(nfloat_ptr res, int sgnbit, gr_ctx_t ctx)

.. function:: int _nfloat_cmp(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
              int _nfloat_cmpabs(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
              int _nfloat_add_1(nfloat_ptr res, ulong x0, slong xexp, int xsgnbit, ulong y0, slong delta, gr_ctx_t ctx)
              int _nfloat_sub_1(nfloat_ptr res, ulong x0, slong xexp, int xsgnbit, ulong y0, slong delta, gr_ctx_t ctx)
              int _nfloat_add_2(nfloat_ptr res, nn_srcptr xd, slong xexp, int xsgnbit, nn_srcptr yd, slong delta, gr_ctx_t ctx)
              int _nfloat_sub_2(nfloat_ptr res, nn_srcptr xd, slong xexp, int xsgnbit, nn_srcptr yd, slong delta, gr_ctx_t ctx)
              int _nfloat_add_3(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
              int _nfloat_sub_3(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
              int _nfloat_add_4(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
              int _nfloat_sub_4(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
              int _nfloat_add_n(nfloat_ptr res, nn_srcptr xd, slong xexp, int xsgnbit, nn_srcptr yd, slong delta, slong nlimbs, gr_ctx_t ctx)
              int _nfloat_sub_n(nfloat_ptr res, nn_srcptr xd, slong xexp, int xsgnbit, nn_srcptr yd, slong delta, slong nlimbs, gr_ctx_t ctx)

Complex numbers
-------------------------------------------------------------------------------

Complex floating-point numbers have the obvious representation as
real pairs.

.. type:: nfloat_complex_ptr
          nfloat_complex_srcptr

.. function:: int nfloat_complex_ctx_init(gr_ctx_t ctx, slong prec, int flags)

.. macro:: NFLOAT_COMPLEX_CTX_DATA_NLIMBS(ctx)
           NFLOAT_COMPLEX_RE(ptr, ctx)
           NFLOAT_COMPLEX_IM(ptr, ctx)
           NFLOAT_COMPLEX_IS_SPECIAL(x, ctx)
           NFLOAT_COMPLEX_IS_ZERO(x, ctx)

.. function:: void nfloat_complex_init(nfloat_complex_ptr res, gr_ctx_t ctx)
              void nfloat_complex_clear(nfloat_complex_ptr res, gr_ctx_t ctx)
              int nfloat_complex_zero(nfloat_complex_ptr res, gr_ctx_t ctx)
              int nfloat_complex_get_acf(acf_t res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_set_acf(nfloat_complex_ptr res, const acf_t x, gr_ctx_t ctx)
              int nfloat_complex_get_acb(acb_t res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_set_acb(nfloat_complex_ptr res, const acb_t x, gr_ctx_t ctx)
              int nfloat_complex_write(gr_stream_t out, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_randtest(nfloat_complex_ptr res, flint_rand_t state, gr_ctx_t ctx)
              void nfloat_complex_swap(nfloat_complex_ptr x, nfloat_complex_ptr y, gr_ctx_t ctx)
              int nfloat_complex_set(nfloat_complex_ptr res, nfloat_complex_ptr x, gr_ctx_t ctx)
              int nfloat_complex_one(nfloat_complex_ptr res, gr_ctx_t ctx)
              int nfloat_complex_neg_one(nfloat_complex_ptr res, gr_ctx_t ctx)
              truth_t nfloat_complex_is_zero(nfloat_complex_srcptr x, gr_ctx_t ctx)
              truth_t nfloat_complex_is_one(nfloat_complex_srcptr x, gr_ctx_t ctx)
              truth_t nfloat_complex_is_neg_one(nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_i(nfloat_complex_ptr res, gr_ctx_t ctx)
              int nfloat_complex_pi(nfloat_complex_ptr res, gr_ctx_t ctx)
              int nfloat_complex_conj(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_re(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_im(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              truth_t nfloat_complex_equal(nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
              int nfloat_complex_set_si(nfloat_complex_ptr res, slong x, gr_ctx_t ctx)
              int nfloat_complex_set_ui(nfloat_complex_ptr res, ulong x, gr_ctx_t ctx)
              int nfloat_complex_set_fmpz(nfloat_complex_ptr res, const fmpz_t x, gr_ctx_t ctx)
              int nfloat_complex_set_fmpq(nfloat_complex_ptr res, const fmpq_t x, gr_ctx_t ctx)
              int nfloat_complex_set_d(nfloat_complex_ptr res, double x, gr_ctx_t ctx)
              int nfloat_complex_set_other(nfloat_complex_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int nfloat_complex_neg(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_add(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
              int nfloat_complex_sub(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
              int _nfloat_complex_sqr_naive(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
              int _nfloat_complex_sqr_standard(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
              int _nfloat_complex_sqr_karatsuba(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
              int _nfloat_complex_sqr(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
              int nfloat_complex_sqr(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int _nfloat_complex_mul_naive(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, nfloat_srcptr c, nfloat_srcptr d, gr_ctx_t ctx)
              int _nfloat_complex_mul_standard(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, nfloat_srcptr c, nfloat_srcptr d, gr_ctx_t ctx)
              int _nfloat_complex_mul_karatsuba(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, nfloat_srcptr c, nfloat_srcptr d, gr_ctx_t ctx)
              int nfloat_complex_mul(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
              int nfloat_complex_mul_2exp_si(nfloat_complex_ptr res, nfloat_complex_srcptr x, slong y, gr_ctx_t ctx)
              int nfloat_complex_cmp(int * res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
              int nfloat_complex_cmpabs(int * res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
              int nfloat_complex_abs(nfloat_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              void _nfloat_complex_vec_init(nfloat_complex_ptr res, slong len, gr_ctx_t ctx)
              void _nfloat_complex_vec_clear(nfloat_complex_ptr res, slong len, gr_ctx_t ctx)
              int _nfloat_complex_vec_zero(nfloat_complex_ptr res, slong len, gr_ctx_t ctx)
              int _nfloat_complex_vec_set(nfloat_complex_ptr res, nfloat_complex_srcptr x, slong len, gr_ctx_t ctx)
              int _nfloat_complex_vec_add(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong len, gr_ctx_t ctx)
              int _nfloat_complex_vec_sub(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong len, gr_ctx_t ctx)
              int nfloat_complex_mat_mul_fixed(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong max_extra_prec, gr_ctx_t ctx)
              int nfloat_complex_mat_mul_block(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong min_block_size, gr_ctx_t ctx)
              int nfloat_complex_mat_mul_reorder(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int nfloat_complex_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int nfloat_complex_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int nfloat_complex_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int nfloat_complex_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)

.. function:: int nfloat_complex_inv(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_div(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
              int nfloat_complex_sqrt(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_rsqrt(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)

    Inverse, division, square root and reciprocal square root (principal
    branches), with errors of a few ulp relative to the magnitude of the
    result; the rounding mode is not respected. Purely real or imaginary
    operands are handled with the real functions.

    The inverse and the quotient are `\bar{y} / |y|^2` and
    `x \bar{y} / |y|^2` with a single reciprocal of `|y|^2` shared by
    both parts (operands whose parts differ greatly in magnitude, or have
    extreme exponents, are rescaled first). At one limb, `|y|^2` is
    formed in registers and truncated to one limb, and both parts of the
    numerator are divided by it with one shared inverse.

    The square root of `z = a + bi` is computed from `r = |z|` (by
    :func:`nfloat_hypot` at one and two limbs and for extreme exponents,
    otherwise as `\sqrt{a^2 + b^2}`, which is cheaper) as `q = \sqrt{(r + |a|)/2}` and `b / (2q)`
    (no cancellation), assigned to the real and imaginary parts according
    to the sign of `a`; the reciprocal square root uses
    `v = 1/\sqrt{2 r^2 (r + |a|)}` and the products `(r + |a|) v` and
    `b v`. At one limb on 64-bit machines, both are computed in registers:
    `r` by the integer square root of the top limbs of `a^2 + b^2`, then
    `1/\sqrt{2(r + |a|)}` (times `1/r` for the reciprocal square root)
    from the one-limb reciprocal square root kernel, and two products.

.. function:: int nfloat_complex_exp(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_expm1(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_exp_pi_i(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_exp2(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_exp10(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_log(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_log1p(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_log_pi_i(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_log2(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_log10(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sin(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_cos(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sin_cos(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_tan(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_cot(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sec(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_csc(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sin_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_cos_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sin_cos_pi(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_tan_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_cot_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sec_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_csc_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sinc(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sinc_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sinh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_cosh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sinh_cosh(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_tanh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_coth(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_sech(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_csch(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_asin(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_acos(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_atan(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_acot(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_asec(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_acsc(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_asinh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_acosh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_atanh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_acoth(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_asech(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_acsch(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_asin_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_acos_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_atan_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_acot_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_asec_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_acsc_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_pow(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)

    Complex elementary functions, as compositions of the real functions.
    Real arguments in the real domain of a function use the real
    function. Otherwise the formulas avoid cancellation except where the
    function itself is ill-conditioned, for example
    `\tan(a+bi) = (\sin a \cos a + i \sinh b \cosh b) /
    (\cos^2 a + \sinh^2 b)` for `|b| < 1/2` and, with `E = \exp(-2|b|)`,
    `\tan(a+bi) = (2E \sin 2a + i \operatorname{sgn}(b) (1 - E^2)) /
    (1 + 2E \cos 2a + E^2)` for `|b| \ge 1/2` (similarly for ``cot``;
    ``sec`` and ``csc`` use the analog of the first formula),
    `\log |z| = \log(a^2 + b^2)/2`, or `\operatorname{log1p}(d)/2` with
    `d = a^2 + b^2 - 1` computed exactly in fixed point when
    `\max(|a|,|b|) \in [1/2, 2)` (also used for `1 - |z|^2` in ``atan``;
    in ``log1p``, `|1+z|^2 - 1 = a(2+a) + b^2` is computed in double
    precision when the terms cancel),
    `\operatorname{atan}(z)` by Kahan-type formulas with ``atan2`` and
    ``log1p``, and ``asin``, ``acos`` and ``acosh`` by Kahan's formulas
    using square roots of `1 \pm z`. The branch cuts are those of
    :type:`acb_t`. The hyperbolic functions and their inverses are
    evaluated as trigonometric functions of `iz`, and ``acot``, ``asec``,
    ``acsc`` (and the hyperbolic versions) as functions of `1/z`.
    There are no internal guard bits: each component typically has an
    error of a few ulp relative to the magnitude of the result (a small
    component of a large result can have a large relative error, and the
    functions of `\pi z` can have an error of order
    `|\pi \operatorname{Im}(z)|` ulp from rounding `\pi \operatorname{Im}(z)`),
    and directed rounding is not supported. Intermediate squares can
    overflow or underflow (giving ``GR_UNABLE``) for arguments near the
    limits of the exponent range. Integer powers `|y| < 2^{16}` are
    computed by binary exponentiation, other powers as `\exp(y \log x)`.

Packed fixed-point arithmetic
-------------------------------------------------------------------------------

A fixed-point number in the range `(-1,1)` with `n`-limb precision
is represented as `n+1` contiguous limbs as follows:

    +---------------+
    |   sign limb   |
    +---------------+
    |  mantissa[0]  |
    +---------------+
    |      ...      |
    +---------------+
    | mantissa[n-1] |
    +---------------+

In the following method signatures, ``nlimbs`` always refers to the
precision ``n`` while the storage is ``nlimbs + 1``.

There is no overflow handling: all methods assume that inputs have
been scaled to a range `[-\varepsilon,\varepsilon]` so that all
intermediate results (including rounding errors) lie in `(-1,1)`.

.. function:: void _nfixed_print(nn_srcptr x, slong nlimbs, slong exp)

    Print the fixed-point number

.. function:: void _nfixed_vec_add(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs)
              void _nfixed_vec_sub(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs)

    Vectorized addition or subtraction of *len* fixed-point numbers.

.. function:: void _nfixed_dot_2(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
              void _nfixed_dot_3(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
              void _nfixed_dot_4(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
              void _nfixed_dot_5(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
              void _nfixed_dot_6(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
              void _nfixed_dot_7(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
              void _nfixed_dot_8(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)

    Dot product with a fixed number of limbs. The ``xstride`` and ``ystride`` parameters
    indicate the offset in number of limbs between consecutive entries
    and may be negative.

.. function:: void _nfixed_mat_mul_classical_precise(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
              void _nfixed_mat_mul_classical(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
              void _nfixed_mat_mul_waksman(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
              void _nfixed_mat_mul_strassen(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong cutoff, slong nlimbs)
              void _nfixed_mat_mul(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)

    Matrix multiplication using various algorithms.
    The *strassen* variant takes a *cutoff* parameter specifying where
    to switch from basecase multiplication to Strassen multiplication.
    The *classical_precise* version computes with one extra limb of
    internal precision; this is only intended for testing purposes.

.. function:: void _nfixed_mat_mul_bound_classical(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong nlimbs)
              void _nfixed_mat_mul_bound_waksman(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong nlimbs)
              void _nfixed_mat_mul_bound_strassen(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong cutoff, slong nlimbs)
              void _nfixed_mat_mul_bound(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong nlimbs)
              void _nfixed_complex_mat_mul_bound(double * bound, double * error, slong m, slong n, slong p, double A, double B, double C, double D, slong nlimbs)

    For the respective matrix multiplication algorithm, computes bounds
    for a size `m \times n \times p` product at precision *nlimbs*
    given entrywise bounds *A* and *B*.

    The *bound* output is set to a bound for the entries in all intermediate
    variables of the computation. This should be < 1 to
    ensure correctness. The *error* output is set to a bound for the
    output error, measured in ulp.
    The caller can assume that the computed bounds are nondecreasing
    functions of *A* and *B*.

    For complex multiplication, the entrywise bounds are for `A+Bi` and `C+Di`.

.. _nfloat-rounding:

Rounding
-------------------------------------------------------------------------------

The rounding mode is set with the context flags ``NFLOAT_RND_FLOOR`` and
``NFLOAT_RND_CEIL`` (see :func:`nfloat_ctx_init`); without either flag, results
are rounded in an unspecified direction for efficiency, usually towards
zero. Warning: directed rounding may not function correctly in combination
with the ``NFLOAT_ALLOW_UNDERFLOW`` flag. In the lists below, the default
mode means neither flag.

**Exact operations.** These return the exact value in every mode:
``set``, ``neg``, ``abs``, ``mul_2exp_si`` (barring over- and underflow),
``floor``, ``ceil``, ``nint``, ``trunc``, ``sgn``, and the integer- or
exact-valued conversions out of ``nfloat`` (``get_arf``, ``get_fmpz``).

**Correctly rounded in every mode.** In the default mode, the result is the
exact value truncated towards zero; with ``NFLOAT_RND_FLOOR`` or
``NFLOAT_RND_CEIL``, the correctly rounded floor or ceiling:

* conversions: ``set_ui``, ``set_si``, ``set_fmpz``, ``set_fmpq``, ``set_d``,
  ``set_arf``, and ``set_other`` with an ``fmpz``, ``fmpq``, ``arf`` or
  ``nfloat`` operand
* ``mul_ui``, ``mul_si``, ``div_ui``, ``div_si``
* ``sqrt``
* ``hypot``, and :func:`nfloat_complex_abs`

**Correctly rounded with directed rounding.** With ``NFLOAT_RND_FLOOR`` or
``NFLOAT_RND_CEIL`` these give the correctly rounded floor or ceiling; in
the default mode the error is as noted:

* ``add``, ``sub``, ``add_ui``, ``add_si``, ``sub_ui``, ``sub_si``: less
  than 2 ulp (the smaller operand is truncated when aligned)
* ``mul``, ``sqr``: the truncation, except that for three or more limbs
  (where an approximate high product is used) it can rarely be one ulp
  smaller in magnitude
* ``div``, ``inv``: the truncation, or for divisors of three or more
  nonzero limbs possibly one ulp larger in magnitude (less than 1 ulp)
* ``rsqrt``: less than `1 + 2^{-37}` ulp
* the vector versions of the above (for example ``_gr_vec_mul_scalar_ui``),
  with the exceptions noted for :func:`_nfloat_vec_div_scalar`

**Valid bounds with directed rounding.** With ``NFLOAT_RND_FLOOR`` or
``NFLOAT_RND_CEIL`` these return a lower respectively upper bound for
the exact value, but not in general the correctly rounded one:

* ``addmul``, ``submul``, ``addmul_ui``, ``addmul_si`` (two roundings, so
  that the result can be far from the exact value, measured in ulp of the
  result, when the terms cancel)
* ``vec_dot``, ``vec_dot_rev``, ``mat_mul``
* the real elementary functions (``exp``, ``log``, ``sin``, ``atan``,
  ``pow``, the less common functions above, etc.): within a few ulp at the
  function precision, and about one ulp in the default mode (see the
  section on elementary functions)

**Not respecting the rounding mode.** These ignore the rounding flags (the
direction of the error is arbitrary):

* fused operations and generic algorithms other than those above
* ``set_str``
* ``set_other`` except for the types listed above
* mixed arithmetic operations with an ``other``, ``fmpz`` or ``fmpq``
  operand
* other transcendental functions (``gamma``, ``zeta``)
* ``nfloat_complex`` operations other than :func:`nfloat_complex_abs`
  (the complex elementary functions have errors of a few ulp relative to
  the magnitude of the result in the default mode)

