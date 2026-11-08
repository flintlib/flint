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
    Powers x^y = exp(y log |x|) by the kernel _mp_real_pow_fp (mp_real.h),
    with the sign (-1)^y for x < 0 and integer y (NaN otherwise); the
    result overflows unless |y log |x|| < 2^(FLINT_BITS - 3).  Simple cases
    (y = +-1, 2, +-1/2, x = 1, powers of two to integer powers) are
    handled directly.
*/

/* 1 if y is an integer, with *odd set to its parity */
static int
_nfloat_is_int(int * odd, nfloat_srcptr y, slong n)
{
    slong e = NFLOAT_EXP(y), b, i;

    if (e <= 0)
        return 0;

    if (e >= FLINT_BITS * n + 1)
    {
        *odd = 0;
        return 1;
    }

    /* the bit of weight 1 is at position b = FLINT_BITS n - e of the mantissa */
    b = FLINT_BITS * n - e;
    for (i = 0; i < b / FLINT_BITS; i++)
        if (NFLOAT_D(y)[i] != 0)
            return 0;
    if (b % FLINT_BITS != 0 && (NFLOAT_D(y)[b / FLINT_BITS] << (FLINT_BITS - b % FLINT_BITS)) != 0)
        return 0;

    *odd = (NFLOAT_D(y)[b / FLINT_BITS] >> (b % FLINT_BITS)) & 1;
    return 1;
}

/* 1 if the mantissa is a power of two */
static int
_nfloat_is_pow2(nfloat_srcptr x, slong n)
{
    return NFLOAT_D(x)[n - 1] == (UWORD(1) << (FLINT_BITS - 1)) && flint_mpn_zero_p(NFLOAT_D(x), n - 1);
}

#if FLINT_BITS == 64
/* the in-register paths of arf/elem.h, rounded into the context */
#define NFLOAT_POW_FAST_DECLINE (-1)

static int
_nfloat_pow_fast_1(nfloat_ptr res, ulong dx, slong ex, ulong dy, slong ey, int ysgn, int rsgn, gr_ctx_t ctx)
{
    ulong y2 = 0, y1 = 0, y0 = 0, err = 0;
    slong k = 0;
    if (_arf_small_pow_1(&y2, &y1, &y0, &k, &err, dx, ex, dy, ey, ysgn) == ARF_ELEM_UNCERTAIN)
        return NFLOAT_POW_FAST_DECLINE;
    return _nfloat_small_set_1(res, y2, y1, y0, k, err, rsgn, ctx);
}

static int
_nfloat_pow_fast_2(nfloat_ptr res, ulong dx1, ulong dx0, slong ex, ulong dy1, ulong dy0,
    slong ey, int ysgn, int rsgn, gr_ctx_t ctx)
{
    ulong y3 = 0, y2 = 0, y1 = 0, y0 = 0, err = 0;
    slong k = 0;
    if (_arf_small_pow_2(&y3, &y2, &y1, &y0, &k, &err, dx1, dx0, ex, dy1, dy0, ey, ysgn) == ARF_ELEM_UNCERTAIN)
        return NFLOAT_POW_FAST_DECLINE;
    return _nfloat_small_set_2(res, y3, y2, y1, y0, k, err, rsgn, ctx);
}
#endif

int
nfloat_pow(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
{
    ulong out[NFLOAT_ELEM_MAX_LIMBS + 2];
    slong n, nk, k;
    ulong err;
    int odd = 0, isint, rsgn;

    if (NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_SPECIAL(y))
    {
        if (NFLOAT_IS_ZERO(y))
            return nfloat_one(res, ctx);
        if (NFLOAT_IS_ZERO(x) && !NFLOAT_IS_SPECIAL(y) && !NFLOAT_SGNBIT(y))
            return nfloat_zero(res, ctx);
        return nfloat_nan(res, ctx);
    }

    n = NFLOAT_CTX_NLIMBS(ctx);
    isint = _nfloat_is_int(&odd, y, n);

    if (NFLOAT_SGNBIT(x) && !isint)
        return nfloat_nan(res, ctx);

    rsgn = NFLOAT_SGNBIT(x) && odd;

    /* |x| = 1 */
    if (NFLOAT_EXP(x) == 1 && _nfloat_is_pow2(x, n))
    {
        nfloat_one(res, ctx);
        NFLOAT_SGNBIT(res) = rsgn;
        return GR_SUCCESS;
    }

    /* y = +-1, 2, +-1/2 */
    if (_nfloat_is_pow2(y, n) && NFLOAT_EXP(y) >= 0 && NFLOAT_EXP(y) <= 2)
    {
        if (NFLOAT_EXP(y) == 1)
            return NFLOAT_SGNBIT(y) ? nfloat_inv(res, x, ctx) : nfloat_set(res, x, ctx);
        if (NFLOAT_EXP(y) == 2 && !NFLOAT_SGNBIT(y))
            return nfloat_sqr(res, x, ctx);
        if (NFLOAT_EXP(y) == 0)
            return NFLOAT_SGNBIT(y) ? nfloat_rsqrt(res, x, ctx) : nfloat_sqrt(res, x, ctx);
    }

    /* (2^e)^y = 2^(e y) for integer y, when it fits: |e|, |y| <
       2^(FLINT_BITS/2 - 2), so that e y + 1 fits in a word */
    if (isint && _nfloat_is_pow2(x, n) && NFLOAT_EXP(y) <= FLINT_BITS / 2 - 2)
    {
        slong e = NFLOAT_EXP(x) - 1;
        slong yi = (slong) (NFLOAT_D(y)[n - 1] >> (FLINT_BITS - NFLOAT_EXP(y)));

        if (NFLOAT_SGNBIT(y))
            yi = -yi;

        if (FLINT_ABS(e) < (WORD(1) << (FLINT_BITS / 2 - 2)))
        {
            nfloat_one(res, ctx);
            NFLOAT_SGNBIT(res) = rsgn;
            NFLOAT_EXP(res) = e * yi + 1;
            NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
            return GR_SUCCESS;
        }
    }

#if FLINT_BITS == 64
    if (n <= 2)
    {
        int status = (n == 1)
            ? _nfloat_pow_fast_1(res, NFLOAT_D(x)[0], NFLOAT_EXP(x), NFLOAT_D(y)[0],
                NFLOAT_EXP(y), NFLOAT_SGNBIT(y), rsgn, ctx)
            : _nfloat_pow_fast_2(res, NFLOAT_D(x)[1], NFLOAT_D(x)[0], NFLOAT_EXP(x),
                NFLOAT_D(y)[1], NFLOAT_D(y)[0], NFLOAT_EXP(y), NFLOAT_SGNBIT(y), rsgn, ctx);
        if (status != NFLOAT_POW_FAST_DECLINE)
            return status;
    }
#endif

    {
        ulong q[3];
        slong qn;
        int tsgn, code;

        nk = _nfloat_elem_limbs(ctx, 0);
        code = _mp_real_pow_fp(out, &err, q, &qn, &tsgn, NFLOAT_D(x), n, NFLOAT_EXP(x),
            NFLOAT_D(y), n, NFLOAT_EXP(y), NFLOAT_SGNBIT(y), _nfloat_elem_bits(ctx, 0), MP_REAL_EXP_FP_MAX_EXP);

        if (code == MP_REAL_POW_HUGE)
            return tsgn ? _nfloat_underflow(res, rsgn, ctx) : _nfloat_overflow(res, rsgn, ctx);
        /* |delta| < 2^-(FLINT_BITS nk + 3), nk from the function
           precision, which may be below the context's */
        if (code == MP_REAL_POW_TINY)
            return _nfloat_one_plus_tiny(res, rsgn, tsgn, -FLINT_BITS * nk - 3, ctx);
        FLINT_ASSERT(code == 0 || code == 1);
        k = (qn == 0) ? 0 : (code ? -(slong) q[0] : (slong) q[0]);
    }

    return _nfloat_set_mpn_err(res, out, nk + 1, k - FLINT_BITS * nk, err, err, rsgn, ctx);
}
