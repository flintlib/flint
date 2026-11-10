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
    Exponentials: the kernels live in arf/mpn_exp.c and the register fast
    paths in arf/elem.h; this file assembles their results into nfloats.
*/

#if FLINT_BITS == 64
static int
_nfloat_exp_fast_1(nfloat_ptr res, ulong d, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong y2, y1, y0;
    slong k = _arf_small_exp_1(&y2, &y1, &y0, d, e, sgnbit);
    return _nfloat_small_set_1(res, y2, y1, y0, k, ARF_SMALL_EXP_ERR_2, 0, ctx);
}

static int
_nfloat_exp_fast_2(nfloat_ptr res, ulong d1, ulong d0, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong y3, y2, y1, y0;
    slong k = _arf_small_exp_2(&y3, &y2, &y1, &y0, d1, d0, e, sgnbit);
    return _nfloat_small_set_2(res, y3, y2, y1, y0, k, ARF_SMALL_EXP_ERR_3, 0, ctx);
}

/* expm1 for -40 <= e <= MP_REAL_EXP_FP_MAX_EXP */
static int
_nfloat_expm1_fast_1(nfloat_ptr res, ulong d, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong y2, y1, y0;
    int sgn;
    slong k = _arf_small_expm1_1(&y2, &y1, &y0, &sgn, d, e, sgnbit);
    return _nfloat_small_set_1(res, y2, y1, y0, k, ARF_SMALL_EXP_ERR_2 + 1, sgn, ctx);
}

static int
_nfloat_expm1_fast_2(nfloat_ptr res, ulong d1, ulong d0, slong e, int sgnbit, gr_ctx_t ctx)
{
    ulong y3, y2, y1, y0;
    int sgn;
    slong k = _arf_small_expm1_2(&y3, &y2, &y1, &y0, &sgn, d1, d0, e, sgnbit);
    return _nfloat_small_set_2(res, y3, y2, y1, y0, k, ARF_SMALL_EXP_ERR_3 + 1, sgn, ctx);
}
#endif

int
nfloat_exp(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong n, e, nk, k;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_pos_inf(res, ctx);
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > MP_REAL_EXP_FP_MAX_EXP)
        return sgnbit ? _nfloat_underflow(res, 0, ctx) : _nfloat_overflow(res, 0, ctx);

    /* |x| < 2^e: exp(x) = 1 + delta, |delta| < 2^(e + 1) <= 2^(-FLINT_BITS n) */
    if (e < -FLINT_BITS * n
        || (e < -NFLOAT_CTX_FUNC_PREC(ctx) && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
        return _nfloat_one_plus_tiny(res, 0, sgnbit, -FLINT_BITS * n, ctx);

#if FLINT_BITS == 64
    if (n == 1)
        return _nfloat_exp_fast_1(res, NFLOAT_D(x)[0], e, sgnbit, ctx);
    if (n == 2)
        return _nfloat_exp_fast_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, sgnbit, ctx);
#endif

    nk = _nfloat_elem_limbs(ctx, 0);
    k = _mp_real_exp_fp(y, &err, NFLOAT_D(x), n, e, sgnbit, _nfloat_elem_bits(ctx, 0));
    return _nfloat_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk, err, err, 0, ctx);
}

int
nfloat_expm1(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    slong n, e, nk, k;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_zero(res, ctx);
        if (NFLOAT_IS_POS_INF(x))
            return nfloat_pos_inf(res, ctx);
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_neg_one(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > MP_REAL_EXP_FP_MAX_EXP)
    {
        if (!sgnbit)
            return _nfloat_overflow(res, 0, ctx);
        /* -1 + exp(x), exp(x) < 2^(-2^(FLINT_BITS - 3)) */
        return _nfloat_one_plus_tiny(res, 1, 1, -FLINT_BITS * n, ctx);
    }

    /* expm1(x) = x (1 + delta), delta = x/2 + x^2/6 + ..., |delta| < |x| */
    if (e < -FLINT_BITS * n
        || (e < -NFLOAT_CTX_FUNC_PREC(ctx) && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
        return _nfloat_set_tiny_rel_err(res, x, sgnbit ? -1 : 1, ctx);

#if FLINT_BITS == 64
    if (n == 1 && e >= -40)
        return _nfloat_expm1_fast_1(res, NFLOAT_D(x)[0], e, sgnbit, ctx);
    if (n == 2 && e >= -40)
        return _nfloat_expm1_fast_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0], e, sgnbit, ctx);
#endif

    {
        ulong elo, ehi;

        /* |x| < 2^-z: an absolute accuracy of 2^-(p + z) */
        nk = _nfloat_elem_limbs(ctx, (e <= 0) ? -e : 0);
        k = _mp_real_expm1_fp(y, &elo, &ehi, NFLOAT_D(x), n, e, sgnbit, _nfloat_elem_bits(ctx, (e <= 0) ? -e : 0));
        return _nfloat_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk, elo, ehi, sgnbit, ctx);
    }
}

int
nfloat_exp2(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong y[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong n, e, nk, k;
    int sgnbit;

    if (NFLOAT_IS_SPECIAL(x))
        return nfloat_exp(res, x, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    sgnbit = NFLOAT_SGNBIT(x);

    if (e > MP_REAL_EXP_FP_MAX_EXP)
        return sgnbit ? _nfloat_underflow(res, 0, ctx) : _nfloat_overflow(res, 0, ctx);

    /* 2^x = 1 + delta, |delta| < 2^e */
    if (e < -FLINT_BITS * n
        || (e < -NFLOAT_CTX_FUNC_PREC(ctx) && !NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)))
        return _nfloat_one_plus_tiny(res, 0, sgnbit, -FLINT_BITS * n, ctx);

    nk = _nfloat_elem_limbs(ctx, 0);
    k = _mp_real_exp2_fp(y, &err, NFLOAT_D(x), n, e, sgnbit, _nfloat_elem_bits(ctx, 0));

    if (err == 0)
    {
        /* an integer: 2^k exactly */
        nfloat_one(res, ctx);
        NFLOAT_EXP(res) = k + 1;
        NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
        return GR_SUCCESS;
    }

    return _nfloat_set_mpn_err(res, y, nk + 1, k - FLINT_BITS * nk, err, err, 0, ctx);
}
