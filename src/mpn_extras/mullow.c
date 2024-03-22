/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*

  0_1_2_3_4_5
0|x x x x x x
1|x x x x x l
2|x x x x l
3|x x x l
4|x x l
5|x l

*/

mp_limb_t
flint_mpn_mullow_basecase(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    mp_limb_t ret;
    mp_size_t ix;

    ret = mpn_mul_1(rp, xp, n, yp[0]);

    for (ix = 1; ix < n; ix++)
    {
        ret += mpn_addmul_1(rp + ix, xp, n - ix, yp[ix]);
        ret += xp[n - ix] * yp[ix];
    }

    return ret;
}
