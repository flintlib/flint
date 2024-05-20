/*
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
void _n_mod_vec_neg(nn_ptr restrict rp, nn_srcptr restrict ip, slong len, ulong mod)
{
    slong ix;

    if (len <= 0)
        FLINT_UNREACHABLE;

    for (ix = 0; ix < len; ix++)
        rp[ix] = (ip[ix] != UWORD(0)) ? mod - ip[ix] : UWORD(0);
}
#if defined(__GNUC__)
# pragma GCC pop_options
#endif
