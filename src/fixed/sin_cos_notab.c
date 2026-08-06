/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fixed.h"

/* fixed_sin_cos_notab: sin(x) and cos(x) for any x in [0, 1)
   without table-based argument reduction.  As in fixed_exp_notab
   the argument is halved h = max(0, r(n) - z) times (z the leading
   zero bits of x, the shift exact within the guard limbs) and
   fixed_sin_cos_reduced runs at depth z + h; the angle is then
   doubled back h times on the cosine alone.  Working with
   g = 1 - cos throughout,

       g(2 theta) = 2 g (2 - g) = 2 (2 g - g^2),

   one sqrhigh per doubling, with every intermediate a pure
   fraction: the doubling chain ends at g(x) <= 1 - cos(1) < 0.46,
   so its inputs stay below g(1/2) < 0.13 and 2g - g^2 below 0.23.
   The final sine comes from one square root,

       sin(x) = sqrt(1 - cos^2) = sqrt(2 g - g^2),

   through mpn_sqrtrem at small sizes and fixed_sqrt_newton above
   the cutoff, the input taken at a limb position of matching parity
   so that the root placement is a limb copy (the same
   normalization-free shape as the burst driver's slice roots).

   Guard limbs: the absolute error of g multiplies by 4 - 4g < 4
   per doubling, and the square root divides by 2 sin(x); with
   x >= 2^(-z-1) on any path that doubles at all, the total
   amplification is bounded by 4^h 2^(z+2) = 2^(2r - z + 2), so
   2h + z + lg n + 8 guard bits cover the doubling chain, the
   sqrhigh slack and the reconstruction together, keeping the
   amplified component below one output ulp.  On the shift-free path
   (h = 0) the outputs are the reduced call's directly.

   The reduction depth r(n): r = 32 measured best across the whole
   rectangular-splitting range on the development machine (r = 24
   and r = 64 trade a few percent inside the noise in spots), with
   r = 24 taking over in the bit-burst regime, matching arb's
   bit-burst sine/cosine (24 halvings).  Retune with the
   forced-depth worker below and profile/p-fixed. */

void
_fixed_sin_cos_notab_r(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int r)
{
    slong xn = n, z, h, G, wn, i, l;
    nn_ptr t, s, g, w, sc;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);

    while (xn > 0 && x[xn - 1] == 0)
        xn--;
    if (xn == 0)
    {
        /* sin(0) = 0, cos(0) = 1 */
        flint_mpn_zero(ysin, n + 1);
        flint_mpn_zero(ycos, n);
        ycos[n] = 1;
        return;
    }

    z = FLINT_BITS * (n - xn)
        + (FLINT_BITS - FLINT_BIT_COUNT(x[xn - 1]));
    h = FLINT_MAX(0, (slong) r - z);

    if (h == 0)
    {
        /* deep enough already: sin and g = 1 - cos directly */
        nn_ptr sg;
        TMP_START;
        sg = TMP_ALLOC((2 * (n + 1)) * sizeof(ulong));
        fixed_sin_cos_reduced(sg, sg + (n + 1), x, n,
            (flint_bitcnt_t) z, 0);
        flint_mpn_copyi(ysin, sg, n);
        ysin[n] = 0;
        ycos[n] = mpn_neg(ycos, sg + (n + 1), n) ? 0 : 1;
        TMP_END;
        return;
    }

    /* guard limbs: 2h + z + lg n + 8 bits (doubling amplification
       4^h, square-root conditioning 2^(z+1), sqrhigh slack), which
       also makes the halving shift exact */
    G = (2 * h + z + FLINT_BIT_COUNT((ulong) n + 128) + 8
        + FLINT_BITS - 1) / FLINT_BITS;
    wn = n + G;

    TMP_START;
    t = TMP_ALLOC((5 * (wn + 2)) * sizeof(ulong));
    s = t + (wn + 2);
    g = s + (wn + 2);
    w = g + (wn + 2);
    sc = w + (wn + 2);

    /* t = (x B^G) 2^-h, exact within the guard */
    {
        slong q = h / FLINT_BITS;
        int b = (int) (h % FLINT_BITS);

        FLINT_ASSERT(q < G || (q == G && b == 0));
        flint_mpn_zero(t, G - q);
        flint_mpn_copyi(t + (G - q), x, n);
        flint_mpn_zero(t + (G - q) + n, q);
        if (b)
            mpn_rshift(t, t, (G - q) + n, b);   /* the outgoing
                bits flow into the guard limbs below x */
    }

    fixed_sin_cos_reduced(s, g, t, wn, (flint_bitcnt_t) (z + h), 0);

    /* double the angle h times on g, then one more 2g - g^2 for the
       sine's square */
    for (i = 0; i <= h; i++)
    {
        flint_mpn_sqrhigh(sc, g, wn);       /* g^2, wn fraction limbs */
        mpn_lshift(w, g, wn, 1);            /* 2g < 1 throughout */
        {
            ulong bw = mpn_sub_n(w, w, sc, wn);
            FLINT_ASSERT(bw == 0);
            (void) bw;
        }
        if (i == h)
            break;                          /* w = sin^2(x) */
        mpn_lshift(g, w, wn, 1);            /* g(2 theta) = 2w < 1 */
    }

    /* cos(x) = 1 - g */
    {
        nn_ptr c = sc;                      /* reuse */
        ulong unit = mpn_neg(c, g, wn) ? 0 : 1;
        flint_mpn_copyi(ycos, c + G, n);
        ycos[n] = unit;
    }

    /* sin(x) = sqrt(w): w in (0, 0.71], with at most about
       2 (r - z) + z + 2 leading zero bits (w ~ x^2 / 2^(2h) scale
       collapses: w = sin(x)^2 >= (x / 2)^2 > 2^(-2z - 4)) */
    l = wn;
    while (l > 1 && w[l - 1] == 0)
        l--;
    if (l == 1 && w[0] == 0)
    {
        flint_mpn_zero(ysin, n + 1);
    }
    else if (wn < FIXED_SIN_COS_NOTAB_SQRT_NEWTON_CUTOFF)
    {
        /* floor(sqrt(w B^(2 wn))): append wn zero limbs; the root
           has ceil((wn + l) / 2) <= wn limbs and is exact */
        nn_ptr V = t, S = s;                /* reuse */
        flint_mpn_zero(V, wn);
        flint_mpn_copyi(V + wn, w, l);
        mpn_sqrtrem(S, NULL, V, wn + l);
        {
            slong sl = (wn + l + 1) / 2;
            flint_mpn_zero(ysin, n + 1);
            FLINT_ASSERT(sl >= G && sl <= wn);
            flint_mpn_copyi(ysin, S + G, sl - G);
        }
    }
    else
    {
        /* Newton at pure limb granularity: take w's top limbs at a
           position vh with vh + wn even, giving vhat in [B^-2, 1),
           fixed_sqrt_newton's accepted range; then
       sin B^wn = sqrt(vhat) B^((vh + wn) / 2) and the placement is
           a limb copy */
        slong nin = wn + 2;
        slong vh = l + ((l + wn) & 1);
        slong off = vh - nin, e, avail;
        nn_ptr vf = t, rt = s;              /* reuse */

        if (off >= 0)
        {
            flint_mpn_copyi(vf, w + off, l - off);
            avail = l - off;
        }
        else
        {
            flint_mpn_zero(vf, -off);
            flint_mpn_copyi(vf - off, w, l);
            avail = l - off;
        }
        if (avail < nin)
            flint_mpn_zero(vf + avail, nin - avail);

        fixed_sqrt_newton(rt, vf, nin, wn + 1);
        /* rt: (wn + 1)-limb fraction of sqrt(vhat) plus a unit limb,
           R as integer; sin B^wn = R B^e with
           e = (vh + wn)/2 - (wn + 1) <= -1 always (sin < 1), so the
           output drops R's low G - e limbs -- a plain copy */
        e = (vh + wn) / 2 - (wn + 1);
        FLINT_ASSERT(e < 0 && G - e <= wn + 2);
        flint_mpn_zero(ysin, n + 1);
        flint_mpn_copyi(ysin, rt + (G - e),
            FLINT_MIN(n + 1, (wn + 2) - (G - e)));
    }

    TMP_END;
}

void
fixed_sin_cos_notab(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n)
{
    _fixed_sin_cos_notab_r(ysin, ycos, x, n,
        (n <= FIXED_SIN_COS_NOTAB_BURST_CUTOFF) ? 32 : 24);
}
