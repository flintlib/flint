/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_mod_vec.h"

#if defined(__GNUC__)
# pragma GCC push_options
# pragma GCC optimize ("unroll-loops")
#endif
void
_n_mod_vec_add_1(nn_ptr restrict rp, nn_srcptr restrict up, nn_srcptr restrict vp, slong len, ulong mod)
{
    slong ix;

    if (len <= 0)
        FLINT_UNREACHABLE;

    for (ix = 0; ix < len; ix++)
    {
        ulong neg = mod - up[ix];
        rp[ix] = (neg > vp[ix]) ? up[ix] + vp[ix] : vp[ix] - neg;
    }
}

void
_n_mod_vec_sub_1(nn_ptr restrict rp, nn_srcptr restrict up, nn_srcptr restrict vp, slong len, ulong mod)
{
    slong ix;

    if (len <= 0)
        FLINT_UNREACHABLE;

    for (ix = 0; ix < len; ix++)
    {
        ulong diff = up[ix] - vp[ix];
        rp[ix] = (up[ix] < vp[ix]) ? mod + diff : diff;
    }
}
#if defined(__GNUC__)
# pragma GCC pop_options
#endif
