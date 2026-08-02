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

/* fixed_exp_notab: exp(x) for any x in [0, 1) without table-based
   argument reduction.  The number of leading zero bits z of x is
   inspected, a target reduction depth r = r(n) is chosen, and the
   argument is halved h = max(0, r - z) times (an exact shift inside
   the guard limbs); fixed_exp_reduced runs on t = x 2^-h at the
   depth z + h, and the result is squared h times,

       exp(x) = exp(x 2^-h)^(2^h),

   each squaring a single sqrhigh at the working precision.  All
   intermediate values exp(x 2^-k) lie in [1, e), so the unit limb
   never carries more than two bits.

   Guard limbs: relative errors double per squaring, so the
   accumulated error is about 2^h (E + h (wn + 3)) ulps of the
   working precision, with E the reduced call's budget and (wn + 3)
   the sqrhigh slack per step; h + lg n + 8 guard BITS cover it, and
   the guard limbs also make the halving shift exact.  With the
   guard in place the amplified component stays below one output
   ulp, leaving the final truncation and, on the shift-free path,
   the reduced call itself as the only visible errors.

   The reduction depth r(n) balances h squarings (~0.4 h M(n))
   against the shorter series of a deeper reduction.  Measured on
   the development machine, r = 32 and r = 16 alternate in windows
   through the rectangular-splitting range -- the r = 16 windows are
   where the automatic dispatch inside fixed_exp_reduced switches to
   the bit-burst path early and beats the series (the "sharp jumps"
   of the series-evaluation cases) -- and r = 32 carries the whole
   bit-burst regime: unlike arb's bit-burst exponential (16
   squarings below ~10^8 bits), the deeper reduction measured
   consistently ahead here, the shrunken first-level splitting trees
   outweighing the extra squarings.  Retune the table with the
   forced-depth worker below and profile/p-fixed. */

static const short _fixed_exp_notab_n_tab[] = { 40, 80, 150, 300 };
static const char _fixed_exp_notab_r_tab[] = { 32, 16, 32, 16, 32 };

void
_fixed_exp_notab_r(nn_ptr y, nn_srcptr x, slong n, int r)
{
    slong xn = n, z, h, G, wn, i;
    nn_ptr t, Y, T;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);

    while (xn > 0 && x[xn - 1] == 0)
        xn--;
    if (xn == 0)
    {
        /* exp(0) = 1 */
        flint_mpn_zero(y, n);
        y[n] = 1;
        return;
    }

    z = FLINT_BITS * (n - xn)
        + (FLINT_BITS - FLINT_BIT_COUNT(x[xn - 1]));
    h = FLINT_MAX(0, (slong) r - z);

    if (h == 0)
    {
        /* already reduced deeply enough: one direct call */
        fixed_exp_reduced(y, x, n, (flint_bitcnt_t) z, 0);
        return;
    }

    /* guard limbs: h + lg n + 8 bits of headroom for the doubling
       relative errors, which also makes the halving shift exact */
    G = (h + FLINT_BIT_COUNT((ulong) n + 128) + 8 + FLINT_BITS - 1)
        / FLINT_BITS;
    wn = n + G;

    TMP_START;
    t = TMP_ALLOC((3 * (wn + 2)) * sizeof(ulong));
    Y = t + (wn + 2);
    T = Y + (wn + 2);

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

    fixed_exp_reduced(Y, t, wn, (flint_bitcnt_t) (z + h), 0);

    /* square h times: Y holds wn fraction limbs and a unit limb; the
       full square spans 2 wn + 2 limbs with a zero top, and sqrhigh
       returns its window's next limb, so the new fixed-point value
       is that returned limb followed by all but the top window limb */
    for (i = 0; i < h; i++)
    {
        T[0] = flint_mpn_sqrhigh(T + 1, Y, wn + 1);
        FLINT_ASSERT(T[wn + 1] <= 1);   /* value below 4 up to slack */
        FLINT_SWAP(nn_ptr, Y, T);
    }

    /* drop the guard limbs */
    flint_mpn_copyi(y, Y + G, n + 1);

    TMP_END;
}

void
fixed_exp_notab(nn_ptr y, nn_srcptr x, slong n)
{
    slong j;
    int r = _fixed_exp_notab_r_tab[
        sizeof(_fixed_exp_notab_n_tab)
            / sizeof(_fixed_exp_notab_n_tab[0])];

    for (j = 0; j < (slong) (sizeof(_fixed_exp_notab_n_tab)
        / sizeof(_fixed_exp_notab_n_tab[0])); j++)
    {
        if (n <= _fixed_exp_notab_n_tab[j])
        {
            r = _fixed_exp_notab_r_tab[j];
            break;
        }
    }
    _fixed_exp_notab_r(y, x, n, r);
}
