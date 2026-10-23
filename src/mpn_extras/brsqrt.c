/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* Up to this many limbs (of the guarded working length) the Newton steps
   use full low products on zero-padded operands. */
#ifndef FLINT_MPN_BRSQRT_MULLOW_CUTOFF
#define FLINT_MPN_BRSQRT_MULLOW_CUTOFF 10
#endif

/*
    2-adic reciprocal square root: develops y with x y^2 == 1 (mod B^n) for
    x == 1 (mod 8), by the Newton iteration

        y' = y - y (x y^2 - 1) / 2.

    This is the p = 2 branch of radix_rsqrtmod_bn with e = FLINT_BITS. Since
    2 is not a unit, the halving is an exact bit shift that loses the top bit
    of the residual, so the iteration is run as a full recompute in bit
    precision with the schedule e_{i+1} = (e_i + 3)/2 (a step with
    x y^2 == 1 (mod 2^c) yields x y'^2 == 1 (mod 2^(2c-2))), and each step
    carries one guard bit above its target precision.

    Entering a step with c bits of equation precision, the low
    zl = (c - 1)/e limbs of x y^2 are 00...01, so its high limbs come from a
    known-low high product; the subtraction of 1 lives entirely below that
    window, the halving shifts in a zero bit from below, and the halved
    residual is divisible by B^zl, so the final product and subtraction act
    only on limbs [zl, nl). In the basecase range everything is done with
    full low products instead.

    The root returned is the one congruent to 1 modulo 4. Returns 1 on
    success, 0 (writing nothing) if x is not 1 modulo 8. res must not alias
    x.
*/

/* keep only the low d bits of the alimbs-limb array a */
static void
_mask(mp_ptr a, mp_size_t d, mp_size_t alimbs)
{
    mp_size_t full = d / FLINT_BITS, r = d % FLINT_BITS;

    if (r)
    {
        a[full] &= (UWORD(1) << r) - 1;
        flint_mpn_zero(a + full + 1, alimbs - full - 1);
    }
    else
    {
        flint_mpn_zero(a + full, alimbs - full);
    }
}

/* single-limb reciprocal square root modulo B; x == 1 (mod 8) */
static mp_limb_t
_n_brsqrt(mp_limb_t x)
{
    mp_limb_t y = 1, t;
    int iter;

    /* each step at least doubles the precision minus two bits:
       3, 4, 6, 10, 18, 34, 66 */
    for (iter = 0; iter < 8; iter++)
    {
        t = x * y * y;
        if (t == 1)
            break;
        t = (t - 1) >> 1;
        y = y - y * t;
    }

    FLINT_ASSERT(x * y * y == 1);
    return y;
}

int
flint_mpn_brsqrt(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)
{
    mp_size_t ed[2 * FLINT_BITS];
    mp_size_t L, k;
    mp_ptr yb, y2, w, t, scr, xp;
    mp_limb_t one = 1;
    TMP_INIT;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(n >= 1);

    if ((x[0] & 7) != 1)
        return 0;

    if (n == 1)
    {
        res[0] = _n_brsqrt(x[0]);
        return 1;
    }

    /* bit-precision schedule from FLINT_BITS n down to the single-limb
       base, which is correct to a full limb */
    ed[0] = FLINT_BITS * n;
    L = 0;
    while (ed[L] > FLINT_BITS)
    {
        ed[L + 1] = (ed[L] + 3) / 2;
        L++;
    }
    ed[L] = FLINT_BITS;

    TMP_START;
    /* all working buffers may need one limb beyond n (the guard bit) */
    yb = TMP_ALLOC(6 * (n + 1) * sizeof(mp_limb_t));
    y2 = yb + (n + 1);
    w = y2 + (n + 1);
    t = w + (n + 1);
    scr = t + (n + 1);
    xp = scr + (n + 1);

    /* x zero-extended to n + 1 limbs, for the padded low products */
    {
        mp_size_t xc = FLINT_MIN(xn, n + 1);
        flint_mpn_copyi(xp, x, xc);
        flint_mpn_zero(xp + xc, n + 1 - xc);
    }

    flint_mpn_zero(yb, n + 1);
    yb[0] = _n_brsqrt(x[0]);

    for (k = L - 1; k >= 0; k--)
    {
        mp_size_t cd = ed[k + 1];              /* x y^2 == 1 (mod 2^cd) */
        mp_size_t nd = ed[k];                  /* target bits */
        mp_size_t hd = nd + 1;                 /* one guard bit */
        mp_size_t nl = (nd + FLINT_BITS - 1) / FLINT_BITS;
        mp_size_t hl = (hd + FLINT_BITS - 1) / FLINT_BITS;
        mp_size_t pl = (cd + FLINT_BITS - 1) / FLINT_BITS;   /* support of y */
        mp_size_t zl = (cd - 1) / FLINT_BITS;  /* zero limbs of the halved residual */

        FLINT_ASSERT(nd <= 2 * cd - 2);
        FLINT_ASSERT(hl <= n + 1);

        if (hl <= FLINT_MPN_BRSQRT_MULLOW_CUTOFF || zl < 1)
        {
            /* y2 = y^2 mod 2^hd */
            flint_mpn_mullow_n(y2, yb, yb, hl);
            _mask(y2, hd, hl);

            /* w = (x y^2 - 1) / 2 mod 2^nd */
            flint_mpn_mullow_n(w, xp, y2, hl);
            _mask(w, hd, hl);
            mpn_sub_1(w, w, hl, 1);
            mpn_rshift(w, w, hl, 1);

            /* y = y - y w mod 2^nd */
            flint_mpn_mullow_n(t, yb, w, nl);
            _mask(t, nd, nl);
            mpn_sub_n(yb, yb, t, nl);
        }
        else
        {
            /* y2 = y^2 mod 2^hd; y has pl nonzero limbs and hl <= 2 pl */
            flint_mpn_mulmid(y2, yb, pl, yb, pl, 0, hl);
            _mask(y2, hd, hl);

            /* w = (x y^2)[zl, hl) = ((x y^2 - 1) / B^zl)[0, hl - zl) */
            _flint_mpn_mulhigh_known_low(w, xp, FLINT_MIN(xn, hl), y2, hl,
                &one, 1, zl, hl, scr);
            _mask(w, hd - zl * FLINT_BITS, hl - zl);
            mpn_rshift(w, w, hl - zl, 1);

            /* y[zl, nl) -= (y w)[0, nl - zl) mod 2^(nd - zl e) */
            flint_mpn_mulmid(t, yb, FLINT_MIN(pl, nl - zl), w, nl - zl, 0, nl - zl);
            _mask(t, nd - zl * FLINT_BITS, nl - zl);
            mpn_sub_n(yb + zl, yb + zl, t, nl - zl);
        }

        _mask(yb, nd, nl);
        if (nl < n + 1)
            flint_mpn_zero(yb + nl, n + 1 - nl);
    }

    flint_mpn_copyi(res, yb, n);

    TMP_END;
    return 1;
}
