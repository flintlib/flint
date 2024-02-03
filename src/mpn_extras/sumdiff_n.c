/*
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

mp_limb_t flint_mpn_sumdiff_n(mp_ptr s, mp_ptr d, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    mp_limb_t ret;
    mp_ptr t;

    if (n == 0)
        return 0;

    if ((s == x && d == y) || (s == y && d == x))
    {
        t = (mp_ptr) flint_malloc(n * sizeof(mp_limb_t));
        ret = mpn_sub_n(t, x, y, n);
        ret += 2 * mpn_add_n(s, x, y, n);
        flint_mpn_copyi(d, t, n);
        flint_free(t);
        return ret;
    }

    if (s == x || s == y)
    {
        ret = mpn_sub_n(d, x, y, n);
        ret += 2 * mpn_add_n(s, x, y, n);
        return ret;
    }

    ret = 2 * mpn_add_n(s, x, y, n);
    ret += mpn_sub_n(d, x, y, n);
    return ret;
}
