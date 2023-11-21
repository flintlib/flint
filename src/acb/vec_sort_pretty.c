/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "acb.h"

#ifndef __compar_fn_t
#if defined(_MSC_VER)
typedef int(*__compar_fn_t) (const void *, const void *);
#else
typedef int(*__compar_fn_t) (__const void *, __const void *);
#endif
#endif

int acb_cmp_pretty(const acb_t a, const acb_t b)
{
    arb_t t, u, v;
    int res;
    arb_init(t);
    arb_init(u);
    arb_init(v);
    arb_abs(u, acb_imagref(a));
    arb_abs(v, acb_imagref(b));
    arb_sub(t, u, v, MAG_BITS);
    res = 0;
    if (arb_contains_zero(t))
    {
        arb_sub(t, acb_realref(a), acb_realref(b), MAG_BITS);
        res = arb_is_positive(t) ? 1 : -1;
    }
    else
    {
        res = arb_is_positive(t) ? 1 : -1;
    }
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    return res;
}

void _acb_vec_sort_pretty(acb_ptr vec, slong len)
{
    qsort(vec, len, sizeof(acb_struct), (__compar_fn_t) acb_cmp_pretty);
}

