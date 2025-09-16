.. _mpn-mod:

**mpn_mod.h** -- integers mod n (packed multi-word n)
===============================================================================

This module provides efficient arithmetic in rings
`R = \mathbb{Z} / n \mathbb{Z}` for medium-sized `n`.
Given an `\ell`-limb modulus `2^{\beta (\ell-1)} \le n < 2^{\beta \ell}`
where `\beta` is ``FLINT_BITS`` (32 or 64),
elements are represented as `\ell`-limb arrays (i.e. ``ulong[l]``),
zero-padded for values that happen to fit in less than `\ell` limbs,
which can be stack-allocated and packed consecutively
without indirection or memory allocation overhead.

This module is designed to use the :ref:`generics <gr>` interface.
As such, the ring is represented by a :type:`gr_ctx_t` context object,
methods return status flags (``GR_SUCCESS``, ``GR_UNABLE``, ``GR_DOMAIN``),
and one can use generic structures such as :type:`gr_poly_t` for
polynomials and :type:`gr_mat_t` for matrices.

Types, macros and constants
-------------------------------------------------------------------------------

.. macro :: MPN_MOD_MIN_LIMBS
            MPN_MOD_MAX_LIMBS

    The number of limbs `\ell` permitted in a modulus. The current limits
    are `2 \le \ell \le 16`, permitting moduli up to 512 bits on
    32-bit machines and 1024 bits on 64-bit machines.
    We exclude single-limb moduli since these are covered by
    :ref:`nmod <nmod>` arithmetic, and this allows not bothering
    with various degenerate cases.
    The upper limit exists so that elements and temporary buffers
    are safe to allocate on the stack and so that simple operations
    like swapping or zeroing elements are not too expensive
    compared to a pointer-and-size representation.
    A second reason is that the algorithms in this module have
    been tuned only for moduli in a certain range.
    For larger moduli, one should use :ref:`fmpz_mod <fmpz-mod>` instead.
    The upper limit might be increased in the future.

Context objects
-------------------------------------------------------------------------------

.. function:: int gr_ctx_init_mpn_mod(gr_ctx_t ctx, const fmpz_t n)
              int _gr_ctx_init_mpn_mod(gr_ctx_t ctx, nn_srcptr n, slong nlimbs)

    Initializes *ctx* to the ring `\mathbb{Z}/n\mathbb{Z}`
    of integers modulo *n* where elements are ``ulong`` arrays with
    the same number of limbs as *n*. This constructor does no
    initialization and returns ``GR_DOMAIN`` if the modulus is nonpositive,
    or ``GR_UNABLE`` if the modulus is not in bounds.

.. function:: void gr_ctx_init_mpn_mod_randtest(gr_ctx_t ctx, flint_rand_t state)

    Initializes *ctx* to a ring with a random modulus.

.. macro:: MPN_MOD_CTX_NLIMBS(ctx)

    Retrieves the number of limbs `\ell` of the modulus.

.. macro:: MPN_MOD_CTX_MODULUS_BITS

    Retrieves the number of bits of the modulus.

.. macro:: MPN_MOD_CTX_MODULUS(ctx)

    Pointer to the limbs of the modulus.

.. macro:: MPN_MOD_CTX_NORM(ctx)

    An integer indicating the number of leading zero bits in the most
    significant limb of the modulus.

.. macro:: MPN_MOD_CTX_MODULUS_NORMED(ctx)

    Pointer to a copy of the modulus left-shifted so that the
    most significant bit is in a limb boundary.

.. macro:: MPN_MOD_CTX_MODULUS_PREINV(ctx)

    Pointer to a precomputed inverse of the (normed) modulus.

.. macro:: MPN_MOD_CTX_IS_PRIME(ctx)

    A :type:`truth_t` flag indicating whether `n` is prime.

.. function:: void mpn_mod_ctx_set_is_field(gr_ctx_t ctx, truth_t is_prime)

    Set the flag indicating whether `n` is prime. Setting this to ``T_TRUE``
    speeds up some algorithms which can assume that the ring
    is actually a field.

Basic operations and arithmetic
-------------------------------------------------------------------------------

.. function:: int mpn_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
              void mpn_mod_ctx_clear(gr_ctx_t ctx)
              truth_t mpn_mod_ctx_is_field(gr_ctx_t ctx)
              void mpn_mod_init(nn_ptr x, gr_ctx_t ctx)
              void mpn_mod_clear(nn_ptr x, gr_ctx_t ctx)
              void mpn_mod_swap(nn_ptr x, nn_ptr y, gr_ctx_t ctx)
              int mpn_mod_set(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
              int mpn_mod_zero(nn_ptr res, gr_ctx_t ctx)
              int mpn_mod_one(nn_ptr res, gr_ctx_t ctx)
              int mpn_mod_set_ui(nn_ptr res, ulong x, gr_ctx_t ctx)
              int mpn_mod_set_si(nn_ptr res, slong x, gr_ctx_t ctx)
              int mpn_mod_neg_one(nn_ptr res, gr_ctx_t ctx)
              int mpn_mod_set_mpn(nn_ptr res, nn_srcptr x, slong xn, gr_ctx_t ctx)
              int mpn_mod_set_fmpz(nn_ptr res, const fmpz_t x, gr_ctx_t ctx)
              int mpn_mod_set_other(nn_ptr res, gr_ptr v, gr_ctx_t v_ctx, gr_ctx_t ctx)
              int mpn_mod_randtest(nn_ptr res, flint_rand_t state, gr_ctx_t ctx)
              int mpn_mod_write(gr_stream_t out, nn_srcptr x, gr_ctx_t ctx)
              int mpn_mod_get_fmpz(fmpz_t res, nn_srcptr x, gr_ctx_t ctx)
              truth_t mpn_mod_is_zero(nn_srcptr x, gr_ctx_t ctx)
              truth_t mpn_mod_is_one(nn_srcptr x, gr_ctx_t ctx)
              truth_t mpn_mod_is_neg_one(gr_srcptr x, gr_ctx_t ctx)
              truth_t mpn_mod_equal(nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
              int mpn_mod_neg(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
              int mpn_mod_add(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
              int mpn_mod_sub(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
              int mpn_mod_add_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx)
              int mpn_mod_sub_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx)
              int mpn_mod_add_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx)
              int mpn_mod_sub_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx)
              int mpn_mod_add_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int mpn_mod_sub_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int mpn_mod_mul(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
              int mpn_mod_mul_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx)
              int mpn_mod_mul_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx)
              int mpn_mod_mul_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int mpn_mod_addmul(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
              int mpn_mod_addmul_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx)
              int mpn_mod_addmul_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx)
              int mpn_mod_addmul_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int mpn_mod_submul(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)
              int mpn_mod_submul_ui(nn_ptr res, nn_srcptr x, ulong y, gr_ctx_t ctx)
              int mpn_mod_submul_si(nn_ptr res, nn_srcptr x, slong y, gr_ctx_t ctx)
              int mpn_mod_submul_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int mpn_mod_fmma(nn_ptr res, nn_srcptr x1, nn_srcptr y1, nn_srcptr x2, nn_srcptr y2, gr_ctx_t ctx)
              int mpn_mod_sqr(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
              int mpn_mod_inv(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
              int mpn_mod_div(nn_ptr res, nn_srcptr x, nn_srcptr y, gr_ctx_t ctx)

    Basic functionality for the ``gr`` method table.
    These methods are interchangeable with their ``gr`` counterparts.
    For example, ``mpn_mod_add(res, x, y, ctx)`` is equivalent to
    ``gr_add(res, x, y, ctx)``.
    The former can be slightly faster as it avoids the indirection of the
    method table lookup.

Vector functions
-------------------------------------------------------------------------------

.. function:: int _mpn_mod_vec_zero(nn_ptr res, slong len, gr_ctx_t ctx)
              int _mpn_mod_vec_clear(nn_ptr res, slong len, gr_ctx_t ctx)
              int _mpn_mod_vec_set(nn_ptr res, nn_srcptr x, slong len, gr_ctx_t ctx)
              void _mpn_mod_vec_swap(nn_ptr vec1, nn_ptr vec2, slong len, gr_ctx_t ctx)
              int _mpn_mod_vec_neg(nn_ptr res, nn_srcptr x, slong len, gr_ctx_t ctx)
              int _mpn_mod_vec_add(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx)
              int _mpn_mod_vec_sub(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx)
              int _mpn_mod_vec_mul(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx)
              int _mpn_mod_vec_mul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx)
              int _mpn_mod_scalar_mul_vec(nn_ptr res, nn_srcptr y, nn_srcptr x, slong len, gr_ctx_t ctx)
              int _mpn_mod_vec_addmul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx)
              int _mpn_mod_vec_submul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx);
              int _mpn_mod_vec_dot(nn_ptr res, nn_srcptr initial, int subtract, nn_srcptr vec1, nn_srcptr vec2, slong len, gr_ctx_t ctx)
              int _mpn_mod_vec_dot_rev(nn_ptr res, nn_srcptr initial, int subtract, nn_srcptr vec1, nn_srcptr vec2, slong len, gr_ctx_t ctx)

    Overrides for generic ``gr`` vector operations with inlined or partially inlined
    code for reduced overhead.

Matrix algorithms
-------------------------------------------------------------------------------

All :type:`gr_mat_t` functionality is supported by this ring.
The following methods implement optimized basic operation overrides
used by higher-level generic routines.

.. function:: int mpn_mod_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

    Waksman's matrix multiplication algorithm using `n^3/2 + O(n^2)` scalar multiplications.
    The operations are done with delayed reduction.

.. function:: int mpn_mod_mat_mul_multi_mod(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

    Reduces matrix multiplication to several ``nmod_mat`` matrix multiplications
    followed by CRT reconstruction. Supports multithreading.

.. function:: int mpn_mod_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

    Dispatches among classical, Waksman and multimodular
    matrix multiplication according to which method is expected
    to perform better for the given dimensions and modulus.
    Strassen is currently not used as the other methods were determined
    to perform better.

.. function:: int mpn_mod_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
              int mpn_mod_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)

    Dispatches to an appropriate generic algorithm (classical
    or block recursive) for triangular solving.

.. function:: int mpn_mod_mat_lu_classical_delayed(slong * res_rank, slong * P, gr_mat_t A, const gr_mat_t A_in, int rank_check, gr_ctx_t ctx)

    Classical LU factorization with delayed modular reductions.

.. function:: int mpn_mod_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)

    Dispatches between classical, delayed-reduction and recursive LU factorization.

.. function:: int mpn_mod_mat_det(nn_ptr res, const gr_mat_t A, gr_ctx_t ctx)

    Dispatches to an appropriate generic algorithm for computing the
    determinant.

Polynomial algorithms
-------------------------------------------------------------------------------

All :type:`gr_poly_t` functionality is supported by this ring.
The following methods implement optimized basic operation overrides
used by higher-level generic routines.

Multiplication
..............

All multiplication algorithms optimize for squaring.

.. function:: int _mpn_mod_poly_mullow_classical(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)

    Polynomial multiplication using the schoolbook algorithm.

.. function:: int _mpn_mod_poly_mullow_KS(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)

    Polynomial multiplication using Kronecker substitution (bit packing).

.. function:: int _mpn_mod_poly_mullow_karatsuba(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, slong cutoff, gr_ctx_t ctx)

    Polynomial multiplication using the Karatsuba algorithm,
    implemented without intermediate modular reductions.
    This algorithm calls itself recursively, switching to
    basecase multiplication (also without intermediate reductions)
    when either *len1* or *len2* is smaller than *cutoff*.

    Currently a full product is computed internally regardless of *len*;
    truncation only skips the modular reductions.

.. function:: int _mpn_mod_poly_mullow_fft_small(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)

    Polynomial multiplication using the small-prime FFT.
    Returns ``GR_UNABLE`` if the small-prime FFT is not available
    or if the coefficients are too large to use this implementation.

.. function:: int _mpn_mod_poly_mullow(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)

    Polynomial multiplication with automatic algorithm selection.

Division
..............

.. function:: int _mpn_mod_poly_inv_series(nn_ptr Q, nn_srcptr B, slong lenB, slong len, gr_ctx_t ctx)
              int _mpn_mod_poly_div_series(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, slong len, gr_ctx_t ctx)

    Power series inversion and division with automatic selection
    between basecase and Newton algorithms.

.. function:: int _mpn_mod_poly_divrem_basecase_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx)
              int _mpn_mod_poly_divrem_basecase(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx)

    Polynomial division with remainder implemented using the basecase
    algorithm with delayed reductions.

.. function:: int _mpn_mod_poly_divrem_q1_preinv1_fmma(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);
              int _mpn_mod_poly_divrem_q1_preinv1_fmma_precond(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);
              int _mpn_mod_poly_divrem_q1_preinv1_karatsuba_precond(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);
              int _mpn_mod_poly_divrem_q1_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);

    Algorithms for polynomial division in the special case where
    `lenA = lenB + 1`. Require `lenB \ge 2`.

.. function:: int _mpn_mod_poly_divrem(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx)
              int _mpn_mod_poly_div(nn_ptr Q, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx)

    Polynomial division with remainder with automatic selection
    between basecase and Newton algorithms.

GCD
..............

.. function:: int _mpn_mod_poly_gcd(nn_ptr G, slong * lenG, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx)

    Polynomial GCD with automatic selection between basecase
    and HGCD algorithms.

.. function:: int _mpn_mod_poly_xgcd(slong * lenG, nn_ptr G, nn_ptr S, nn_ptr T, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, gr_ctx_t ctx);

    Polynomial extended GCD with automatic selection between basecase
    and HGCD algorithms.
