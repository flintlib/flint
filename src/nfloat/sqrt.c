/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "gr.h"
#include "nfloat.h"
#include "arf/mpn_small.h"
#include "impl.h"

/*
    Square root and reciprocal square root on the root kernels of
    arf (arf/mpn_div.c): for x = A 2^e with A = a / B^n in [1/2, 1) and t = e mod 2,
    sqrt(x) = sqrt(A / 2^t) 2^((e + t) / 2) and 1/sqrt(x) =
    M 2^((1 - t) - (e - t) / 2) with M = (2^(2-t) A)^(-1/2), the kernels
    returning the mantissas in [1/2, 1) with the corrections 0 or 1 to
    these exponents. Without directed rounding the _approx kernels (the
    truncated root; the reciprocal root within 1 + 2^-37 ulp), with
    directed rounding the rounded kernels toward or away from zero.
*/

FLINT_FORCE_INLINE int
_nfloat_root_rnd(gr_ctx_t ctx)
{
    return nfloat_should_round_up(0, ctx) ? ARF_RND_UP : ARF_RND_DOWN;
}

int
nfloat_sqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n, e, fix;
    int t;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_nan(res, ctx);
        else
            return nfloat_set(res, x, ctx);
    }

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    t = e & 1;

    if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
    {
#if FLINT_BITS == 64
        if (n == 1)
            _arf_mpn_sqrt_1_approx(NFLOAT_D(res), &fix, NFLOAT_D(x)[0], t);
        else if (n == 2)
            _arf_mpn_sqrt_2_approx(NFLOAT_D(res), &fix, NFLOAT_D(x)[1], NFLOAT_D(x)[0], t);
        else
#endif
        {
            /* the exact truncated root of the mantissa padded to 2n words
               (what the _approx kernel does below the Newton range, which
               lies beyond NFLOAT_MAX_LIMBS: called directly) */
            flint_mpn_sqrthigh_exact(NFLOAT_D(res), NFLOAT_D(x), n, t, n - 1);
            fix = 0;
        }
    }
    else
    {
        int rnd = _nfloat_root_rnd(ctx);

#if FLINT_BITS == 64
        if (n == 1)
            _arf_mpn_sqrt_1(NFLOAT_D(res), &fix, NFLOAT_D(x)[0], t, rnd);
        else if (n == 2)
            _arf_mpn_sqrt_2(NFLOAT_D(res), &fix, NFLOAT_D(x)[1], NFLOAT_D(x)[0], t, rnd);
        else
#endif
            _arf_mpn_sqrt(NFLOAT_D(res), &fix, NFLOAT_D(x), n, t, n * FLINT_BITS, rnd);
    }

    NFLOAT_EXP(res) = (e + t) / 2 + fix;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}

int
nfloat_rsqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n, e, fix;
    nn_srcptr a;
    int t;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_pos_inf(res, ctx);
        else if (NFLOAT_IS_POS_INF(x))
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx);
    }

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    t = e & 1;
    a = NFLOAT_D(x);

    /* x = 2^(e - 1) with e - 1 even (the small kernels and the Newton
       iteration require M < 1) */
    if (t && a[n - 1] == (UWORD(1) << (FLINT_BITS - 1)) && flint_mpn_zero_p(a, n - 1))
    {
        if (res != x)
            _nfloat_copy_limbs(NFLOAT_D(res), a, n);
        NFLOAT_EXP(res) = (1 - e) / 2 + 1;
        NFLOAT_SGNBIT(res) = 0;
        return GR_SUCCESS;
    }

    if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
    {
#if FLINT_BITS == 64
        if (n == 1)
            _arf_mpn_rsqrt_1_approx(NFLOAT_D(res), &fix, a[0], t);
        else if (n == 2)
            _arf_mpn_rsqrt_2_approx(NFLOAT_D(res), &fix, a[1], a[0], t);
        else
#endif
        if (n + 1 <= FLINT_MPN_RSQRTHIGH_MP_REAL_CUTOFF)
        {
            /* the Newton iteration with a guard word (what the _approx
               kernel does in this range: called directly) */
            ulong Y[NFLOAT_MAX_LIMBS + 1];
            _flint_mpn_rsqrthigh_newton(Y, n + 1, a, n, t);
            if (FLINT_UNLIKELY(!(Y[n] >> (FLINT_BITS - 1))))
            {
                /* the approximation dipped below 1/2: clamped */
                _nfloat_zero_limbs(Y + 1, n - 1);
                Y[n] = UWORD(1) << (FLINT_BITS - 1);
            }
            _nfloat_copy_limbs(NFLOAT_D(res), Y + 1, n);
            fix = 0;
        }
        else
            _arf_mpn_rsqrt_approx(NFLOAT_D(res), &fix, a, n, t, n * FLINT_BITS);
    }
    else
    {
        int rnd = _nfloat_root_rnd(ctx);

#if FLINT_BITS == 64
        if (n == 1)
            _arf_mpn_rsqrt_1(NFLOAT_D(res), &fix, a[0], t, rnd);
        else if (n == 2)
            _arf_mpn_rsqrt_2(NFLOAT_D(res), &fix, a[1], a[0], t, rnd);
        else
#endif
            _arf_mpn_rsqrt(NFLOAT_D(res), &fix, a, n, t, n * FLINT_BITS, rnd);
    }

    NFLOAT_EXP(res) = (1 - t) - (e - t) / 2 + fix;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}
