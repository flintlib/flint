/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    2-adic square root modulo B^n by one Karp-Markstein refinement on top of
    a half-precision reciprocal square root; the p = 2 branch of
    radix_sqrtmod_bn with e = FLINT_BITS. With m = ceil(n/2) and mm = m + 1:

        y = x^(-1/2) mod B^mm             (reciprocal square root)
        b = x y mod B^mm                  (b^2 == x mod B^mm)
        s = b + y (x - b^2) / 2 mod B^n   (one Karp-Markstein step)

    Since x - b^2 == 0 (mod 2^(e mm)) is even, the halving is exact but is
    only guaranteed to be aligned at bit e mm - 1, so the correction

        B^(mm-1) ((y dh) 2^(e-1) mod B^(n-mm+1)),   dh = (x - b^2)[mm, n+1),

    straddles a limb boundary and is added over limbs [mm-1, n) of the root.
    The reciprocal square root is taken one limb beyond ceil(n/2) so that the
    two bits lost by the 2-adic squaring leave all n limbs correct.

    Returns the root congruent to 1 modulo 4 (1 on success), or 0 if x is
    not 1 modulo 8. res must not alias x.
*/
int
flint_mpn_bsqrt(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)
{
    mp_size_t m, nm, mm;
    mp_ptr y, b;
    TMP_INIT;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(n >= 1);

    if ((x[0] & 7) != 1)
        return 0;

    m = (n + 1) / 2;
    nm = n - m;
    mm = m + 1;

    TMP_START;
    y = TMP_ALLOC((2 * mm + 2 * nm + (n + 1)) * sizeof(mp_limb_t));
    b = y + mm;

    flint_mpn_brsqrt(y, x, xn, mm);

    /* b = x y mod B^mm */
    flint_mpn_mulmid(b, x, FLINT_MIN(xn, mm), y, mm, 0, mm);
    flint_mpn_copyi(res, b, FLINT_MIN(mm, n));

    /* nm = n - m = (n + 1) - mm is both the length of (x - b^2)[mm, n+1)
       and the length of the correction window [mm-1, n) */
    if (nm > 0)
    {
        mp_size_t w = n + 1, ah;
        mp_ptr dh = b + mm, scr = dh + nm, corr = scr + w;

        /* dh = (b^2)[mm, n+1); the low mm limbs of b^2 are x mod B^mm */
        _flint_mpn_mulhigh_known_low(dh, b, mm, b, mm, x, FLINT_MIN(xn, mm),
            mm, w, scr);

        /* dh = (x - b^2)[mm, n+1), x being zero above limb xn */
        mpn_neg(dh, dh, nm);
        ah = (xn > mm) ? FLINT_MIN(xn - mm, nm) : 0;
        if (ah > 0)
            mpn_add(dh, dh, nm, x + mm, ah);

        /* s[mm-1, n) += (y dh 2^(e-1))[0, nm); only the low nm <= mm limbs
           of y enter, and the carry out is the reduction mod B^n */
        flint_mpn_mulmid(corr, y, FLINT_MIN(mm, nm), dh, nm, 0, nm);
        mpn_lshift(corr, corr, nm, FLINT_BITS - 1);

        if (n > mm)
            flint_mpn_zero(res + mm, n - mm);
        mpn_add_n(res + (mm - 1), res + (mm - 1), corr, nm);
    }

    TMP_END;
    return 1;
}
