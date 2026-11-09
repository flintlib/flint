/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "elem.h"

/*
    Hyperbolic functions: the kernels _mp_real_sinh_cosh_fp and
    _mp_real_tanh_fp (mp_real.h) on the mantissa, after the tiny
    arguments and the in-register fast paths of arf/elem.h at one and two
    limbs; coth, sech and csch as reciprocals (_nfloat_set_mpn_err_recip).
*/

/* the buffers of the kernels: w <= nk + (-e + 1) / FLINT_BITS + 2 for
   |x| >= 2^(-FLINT_BITS n / 2) (one more limb for coth) */
#define HYP_LIMBS (NFLOAT_ELEM_MAX_LIMBS + 2)

/* sinh (rs), cosh (rc), either NULL */
static int
_nfloat_sinh_cosh(nfloat_ptr rs, nfloat_ptr rc, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong ys[HYP_LIMBS], yc[HYP_LIMBS];
    slong n, e, w, k;
    ulong err;
    int sgnbit, status = GR_SUCCESS;

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > MP_REAL_EXP_FP_MAX_EXP)
    {
        if (rs != NULL)
            status |= _nfloat_overflow(rs, sgnbit, ctx);
        if (rc != NULL)
            status |= _nfloat_overflow(rc, 0, ctx);
        return status;
    }

    /* sinh x = x (1 + delta), cosh x = 1 + delta', delta, delta' < x^2 */
    if (2 * e <= -FLINT_BITS * n)
    {
        if (rs != NULL)
            status |= _nfloat_set_tiny_rel_err(rs, x, 1, ctx);
        if (rc != NULL)
            status |= _nfloat_one_plus_tiny(rc, 0, 0, -FLINT_BITS * n, ctx);
        return status;
    }

#if FLINT_BITS == 64
    if (n == 1 && (rs == NULL || e >= -31))
    {
        ulong s2 = 0, s1 = 0, s0 = 0, c2 = 0, c1 = 0, c0 = 0;
        _arf_small_sinh_cosh_1((rs != NULL) ? &s2 : NULL, &s1, &s0, (rc != NULL) ? &c2 : NULL, &c1, &c0,
            &k, NFLOAT_D(x)[0], e);
        if (rs != NULL)
            status |= _nfloat_small_set_1(rs, s2, s1, s0, k, ARF_SMALL_HYP_ERR_2, sgnbit, ctx);
        if (rc != NULL)
            status |= _nfloat_small_set_1(rc, c2, c1, c0, k, ARF_SMALL_HYP_ERR_2, 0, ctx);
        return status;
    }
    if (n == 2 && (rs == NULL || e >= -40))
    {
        ulong s3 = 0, s2 = 0, s1 = 0, s0 = 0, c3 = 0, c2 = 0, c1 = 0, c0 = 0;
        _arf_small_sinh_cosh_2((rs != NULL) ? &s3 : NULL, &s2, &s1, &s0, (rc != NULL) ? &c3 : NULL, &c2, &c1, &c0,
            &k, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e);
        if (rs != NULL)
            status |= _nfloat_small_set_2(rs, s3, s2, s1, s0, k, ARF_SMALL_HYP_ERR_3, sgnbit, ctx);
        if (rc != NULL)
            status |= _nfloat_small_set_2(rc, c3, c2, c1, c0, k, ARF_SMALL_HYP_ERR_3, 0, ctx);
        return status;
    }
#endif

    k = _mp_real_sinh_cosh_fp((rs != NULL) ? ys : NULL, (rc != NULL) ? yc : NULL, &w, &err,
        NFLOAT_D(x), n, e, _nfloat_elem_bits(ctx, 0));

    if (rs != NULL)
        status |= _nfloat_set_mpn_err(rs, ys, w + 1, k - FLINT_BITS * w, err, err, sgnbit, ctx);
    if (rc != NULL)
        status |= _nfloat_set_mpn_err(rc, yc, w + 1, k - FLINT_BITS * w, err, err, 0, ctx);
    return status;
}

int
nfloat_sinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_set(res, x, ctx) | (NFLOAT_IS_NAN(x) ? nfloat_nan(res, ctx) : GR_SUCCESS);

    return _nfloat_sinh_cosh(res, NULL, x, ctx);
}

int
nfloat_cosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_INF(x))
            return nfloat_pos_inf(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_sinh_cosh(NULL, res, x, ctx);
}

int
nfloat_sinh_cosh(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        /* sinh first (it leaves special values unchanged, so res1 may
           alias x), as two statements: the operands of | are
           unsequenced */
        int status = nfloat_sinh(res1, x, ctx);
        return status | nfloat_cosh(res2, x, ctx);
    }

    if (res1 == x || res2 == x)
    {
        ulong t[NFLOAT_MAX_ALLOC];
        nfloat_set(t, x, ctx);
        return _nfloat_sinh_cosh(res1, res2, t, ctx);
    }

    return _nfloat_sinh_cosh(res1, res2, x, ctx);
}

/* tanh, coth, sech, csch (which = ARF_SMALL_HYP_*) of a normal x by
   the register paths (2^-32 <= |x| < 2^61 resp. 2^-41), or
   NFLOAT_FAST_DECLINE */
static int
_nfloat_hyp_fast(nfloat_ptr res, nfloat_srcptr x, int which, gr_ctx_t ctx)
{
#if FLINT_BITS == 64
    slong n = NFLOAT_CTX_NLIMBS(ctx), e = NFLOAT_EXP(x), k;
    ulong err;
    int sgnbit = (which == ARF_SMALL_HYP_SECH) ? 0 : NFLOAT_SGNBIT(x);

    if (n > 2 || e > MP_REAL_EXP_FP_MAX_EXP || e < ((n == 1) ? -31 : -40))
        return NFLOAT_FAST_DECLINE;

    if (n == 1)
    {
        ulong t2, t1, t0;
        if (_arf_small_tanh_1(&t2, &t1, &t0, &k, &err, NFLOAT_D(x)[0], e, which) == ARF_SMALL_HYP_ONE)
            return _nfloat_small_set_1(res, UWORD(1), UWORD(0), UWORD(0), 0, 1, sgnbit, ctx);
        return _nfloat_small_set_1(res, t2, t1, t0, k, err, sgnbit, ctx);
    }
    else
    {
        ulong t3, t2, t1, t0;
        if (_arf_small_tanh_2(&t3, &t2, &t1, &t0, &k, &err, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, which) == ARF_SMALL_HYP_ONE)
            return _nfloat_small_set_2(res, UWORD(1), UWORD(0), UWORD(0), UWORD(0), 0, 1, sgnbit, ctx);
        return _nfloat_small_set_2(res, t3, t2, t1, t0, k, err, sgnbit, ctx);
    }
#else
    return NFLOAT_FAST_DECLINE;
#endif
}

/* tanh (recip = 0) or coth (recip = 1) of a normal x */
static int
_nfloat_tanh_coth(nfloat_ptr res, nfloat_srcptr x, int recip, gr_ctx_t ctx)
{
    ulong q[HYP_LIMBS];
    slong n, e, w;
    ulong err;
    int sgnbit, status;

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    /* tanh x = x (1 - delta), coth x = (1/x) (1 + delta), delta < x^2 / 3
       (for coth three more bits: 1/x is itself within an ulp, so delta
       must stay below an eighth of one for the result to stay within
       about an ulp) */
    if (2 * e <= -FLINT_BITS * n - (recip ? 3 : 0))
    {
        if (!recip)
            return _nfloat_set_tiny_rel_err(res, x, -1, ctx);
        /* 1/x within an ulp, then the relative error delta < 2^(2e - 1) */
        status = nfloat_inv(res, x, ctx);
        if (status == GR_SUCCESS && !NFLOAT_IS_SPECIAL(res))
            status = _nfloat_set_tiny_rel_err(res, res, 1, ctx);
        return status;
    }

    status = _nfloat_hyp_fast(res, x, recip ? ARF_SMALL_HYP_COTH : ARF_SMALL_HYP_TANH, ctx);
    if (status != NFLOAT_FAST_DECLINE)
        return status;

    /* |tanh x| = 1 - 2 / (E^2 + 1) with 2 / (E^2 + 1) < 2 exp(-2 |x|) <
       2^(1 - 2.88 |x|) <= 2^(-FLINT_BITS n) for |x| >= 2^(e - 1) >= 23 n + 1
       (and |coth x| = 1 + 2 / (E^2 - 1) likewise) */
    if (e - 1 >= (slong) FLINT_BIT_COUNT(23 * n + 1))
        return _nfloat_one_plus_tiny(res, sgnbit, !recip, -FLINT_BITS * n, ctx);

    /* (coth: one more limb, so that the quotient has its wn limbs
       without scaling up the error: _nfloat_set_mpn_err_recip) */
    _mp_real_tanh_fp(q, &w, &err, NFLOAT_D(x), n, e, _nfloat_elem_bits(ctx, 0) + (recip ? FLINT_BITS : 0));
    return _nfloat_set_mpn_err_recip(res, q, w + 1, -FLINT_BITS * w, err, recip, sgnbit, ctx);
}

int
nfloat_tanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_neg_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_tanh_coth(res, x, 0, ctx);
}

int
nfloat_coth(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_neg_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_tanh_coth(res, x, 1, ctx);
}

/* sech (which = ARF_SMALL_HYP_SECH) or csch of a normal x: the
   reciprocal of cosh resp. sinh */
static int
_nfloat_sech_csch(nfloat_ptr res, nfloat_srcptr x, int which, gr_ctx_t ctx)
{
    ulong y[HYP_LIMBS];
    slong n, e, w, k;
    ulong err;
    int sgnbit, status, csch = (which == ARF_SMALL_HYP_CSCH);

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = csch ? NFLOAT_SGNBIT(x) : 0;

    /* sech x, |csch x| < 2 exp(-|x|) underflow for |x| >= 2^(FLINT_BITS - 3) */
    if (e > MP_REAL_EXP_FP_MAX_EXP)
        return _nfloat_underflow(res, sgnbit, ctx);

    /* sech x = 1 - delta, 0 < delta < x^2 / 2; csch x = (1/x)(1 - delta),
       0 < delta < x^2 / 6 */
    if (2 * e <= -FLINT_BITS * n)
    {
        if (!csch)
            return _nfloat_one_plus_tiny(res, 0, 1, -FLINT_BITS * n, ctx);
        status = nfloat_inv(res, x, ctx);
        if (status == GR_SUCCESS && !NFLOAT_IS_SPECIAL(res))
            status = _nfloat_set_tiny_rel_err(res, res, -1, ctx);
        return status;
    }

    status = _nfloat_hyp_fast(res, x, which, ctx);
    if (status != NFLOAT_FAST_DECLINE)
        return status;

    /* (csch: one more limb, as for coth) */
    k = _mp_real_sinh_cosh_fp(csch ? y : NULL, csch ? NULL : y, &w, &err,
        NFLOAT_D(x), n, e, _nfloat_elem_bits(ctx, 0) + (csch ? FLINT_BITS : 0));
    return _nfloat_set_mpn_err_recip(res, y, w + 1, k - FLINT_BITS * w, err, 1, sgnbit, ctx);
}

int
nfloat_sech(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_INF(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_sech_csch(res, x, ARF_SMALL_HYP_SECH, ctx);
}

int
nfloat_csch(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_INF(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    return _nfloat_sech_csch(res, x, ARF_SMALL_HYP_CSCH, ctx);
}
