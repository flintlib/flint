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
    Windowed middle product via a balanced low product.  The low zhi limbs of
    a * b are (a mod B^zhi) * (b mod B^zhi) mod B^zhi, computed exactly by
    flint_mpn_mullow_n on operands zero-padded / truncated to n = zhi limbs; the
    window [zlo, zhi) is then the top zn of those.  Correct for any input; the
    result is the *exact* window (with carry-in from below zlo), a zero-deficit
    instance of the contract.  Cheapest when zhi is small (a low or nearly-low
    window); for large zhi it computes an almost full n x n product.
*/
void
flint_mpn_mulmid_via_mullow_n(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn,
                              mp_size_t zlo, mp_size_t zhi)
{
    mp_size_t n = zhi;
    mp_size_t zn = zhi - zlo;
    mp_srcptr xp = a, yp = b;
    mp_ptr rp;
    TMP_INIT;

    FLINT_ASSERT(an >= 1 && bn >= 1);
    FLINT_ASSERT(0 <= zlo && zlo < zhi && zhi <= an + bn);

    TMP_START;

    /* present both operands as n limbs: truncate if longer, zero-extend if shorter */
    if (an < n)
    {
        mp_ptr X = TMP_ARRAY_ALLOC(n, mp_limb_t);
        flint_mpn_copyi(X, a, an);
        flint_mpn_zero(X + an, n - an);
        xp = X;
    }
    if (bn < n)
    {
        mp_ptr Y = TMP_ARRAY_ALLOC(n, mp_limb_t);
        flint_mpn_copyi(Y, b, bn);
        flint_mpn_zero(Y + bn, n - bn);
        yp = Y;
    }

    if (zlo == 0)                       /* window is exactly the low n limbs */
    {
        flint_mpn_mullow_n(z, xp, yp, n);
        TMP_END;
        return;
    }

    rp = TMP_ARRAY_ALLOC(n, mp_limb_t);
    flint_mpn_mullow_n(rp, xp, yp, n);
    flint_mpn_copyi(z, rp + zlo, zn);
    TMP_END;
}
