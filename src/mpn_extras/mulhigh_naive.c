/*
    Copyright (C) 2024 Albin Ahlbäck
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

mp_limb_t
_flint_mpn_mulhigh_n_naive(mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n)
{
    mp_limb_t low;
    mp_limb_t a, b;

    if (n == 1)
    {
        umul_ppmm(rp[0], low, up[0], vp[0]);
    }
    else if (n == 2)
    {
        FLINT_MPN_MUL_2X2(rp[1], rp[0], low, b, up[1], up[0], vp[1], vp[0]);
    }
    else if (n == 3)
    {
        NN_DOTREV_S3_1X1_HIGH(b, a, up, vp, 2);
        NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, up, vp, 3);
        NN_DOTREV_S3_A3_1X1(b, a, rp[0], 0, b, a, up + 1, vp + 1, 2);
        NN_ADDMUL_S2_A2_1X1(rp[2], rp[1], b, a, up[2], vp[2]);
    }
    else
    {
        mp_limb_t tp[2];
        slong ix;

        NN_DOTREV_S3_1X1_HIGH(b, a, up, vp, n - 1);
        NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, up, vp, n);
        tp[0] = a;
        tp[1] = b;

        umul_ppmm(rp[1], rp[0], up[n - 1], vp[1]);
        for (ix = 2; ix < n; ix++)
            rp[ix] = mpn_addmul_1(rp, up + n - ix, ix, vp[ix]);

        mpn_add(rp, rp, n, tp, 2);
    }

    return low;
}
