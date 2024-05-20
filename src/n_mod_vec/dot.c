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

#include "n_mod.h"
#include "n_mod_vec.h"

#if defined(__GNUC__)
# pragma GCC push_options
# pragma GCC optimize ("unroll-loops")
#endif
ulong
_n_mod_vec_dot_0(nn_srcptr up, nn_srcptr vp, slong len, n_mod_ctx_srcptr ctx)
{
    ulong r0 = 0;
    slong ix;

    if (len <= 0)
        FLINT_UNREACHABLE;

    for (ix = 0; ix < len; ix++)
        r0 += up[ix] * vp[ix];

    return n_mod_set_ui(r0, ctx);
}

ulong
_n_mod_vec_dot_1(nn_srcptr up, nn_srcptr vp, slong len, n_mod_ctx_srcptr ctx)
{
    ulong r0 = 0, r1 = 0;
    slong ix;

    if (len <= 0)
        FLINT_UNREACHABLE;

    for (ix = 0; ix < len; ix++)
    {
        ulong prod = up[ix] * vp[ix];

#if 0
        ulong r0_old = r0;
        r0 += prod;
        r1 += (r0 < r0_old);
#elif defined(__GNUC__)
        r1 += __builtin_add_overflow(prod, r0, &r0);
#else
        add_ssaaaa(r1, r0, r1, r0, 0, prod);
#endif
    }

    return n_mod_set_uiui(r1, r0, ctx);
}
#if defined(__GNUC__)
# pragma GCC pop_options
#endif
