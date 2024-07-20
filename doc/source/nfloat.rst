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
Atomic operations typically give slightly worse than correct rounding,
e.g. with 1-2 ulp error. The rounding is not guaranteed to be identical
on 64-bit and 32-bit machines.
Planned features include:

* Directed rounding modes
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

    Returns ``GR_UNABLE`` without initializating the context object
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

Infinities and NaNs are disabled by default to improve performance,
as this allows certain functions to skip checks for such values.

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

.. function:: int nfloat_sqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_rsqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

.. function:: int nfloat_sgn(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_im(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

.. function:: int nfloat_floor(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_ceil(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_trunc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_nint(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

.. function:: int nfloat_pow(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

.. function:: int nfloat_pi(nfloat_ptr res, gr_ctx_t ctx)
              int nfloat_exp(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_expm1(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_log(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_log1p(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_cos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_tan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_sinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_cosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_tanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_atan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_gamma(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
              int nfloat_zeta(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

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

.. function:: int _nfloat_vec_dot(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
              int _nfloat_vec_dot_rev(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)

Matrix functions
-------------------------------------------------------------------------------

.. function:: int nfloat_mat_mul_fixed_classical(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int nfloat_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int nfloat_mat_mul_block(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong min_block_size, gr_ctx_t ctx)
              int nfloat_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

    Different implementations of matrix multiplication.

.. function:: int nfloat_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int nfloat_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int nfloat_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)

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
              int nfloat_complex_inv(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
              int nfloat_complex_div(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
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
              int nfloat_complex_mat_mul_fixed_classical(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int nfloat_complex_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int nfloat_complex_mat_mul_block(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong min_block_size, gr_ctx_t ctx)
              int nfloat_complex_mat_mul_reorder(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int nfloat_complex_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int nfloat_complex_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int nfloat_complex_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int nfloat_complex_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
