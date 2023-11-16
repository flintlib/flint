/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

/* this can be improved */
static int
use_recurrence(const acb_t n, const acb_t m, slong prec)
{
    if (!acb_is_int(n) || !arb_is_nonnegative(acb_realref(n)))
        return 0;

    if (arf_cmpabs_ui(arb_midref(acb_realref(n)), prec) > 0)
        return 0;

    if (arf_cmpabs(arb_midref(acb_realref(n)), arb_midref(acb_realref(m))) >= 0)
        return 0;

    return 1;
}

void
acb_hypgeom_laguerre_l_ui_recurrence(acb_t res, ulong n, const acb_t m,
    const acb_t z, slong prec)
{
    acb_t t, u, v;
    ulong k;

    if (n == 0)
    {
        acb_one(res);
        return;
    }

    if (n == 1)
    {
        acb_sub(res, m, z, prec);
        acb_add_ui(res, res, 1, prec);
        return;
    }

    acb_init(t);
    acb_init(u);
    acb_init(v);

    acb_one(t);
    acb_sub(u, m, z, prec);
    acb_add_ui(u, u, 1, prec);

    for (k = 2; k <= n; k++)
    {
        acb_add_ui(v, m, k - 1, prec);
        acb_mul(t, t, v, prec);

        acb_add_ui(v, m, 2 * k - 1, prec);
        acb_sub(v, v, z, prec);
        acb_mul(v, v, u, prec);

        acb_sub(t, v, t, prec);
        acb_div_ui(t, t, k, prec);

        acb_swap(t, u);
    }

    acb_set(res, u);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
acb_hypgeom_laguerre_l(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec)
{
    acb_t t, u, v;

    if (use_recurrence(n, m, prec))
    {
        acb_hypgeom_laguerre_l_ui_recurrence(res,
            arf_get_si(arb_midref(acb_realref(n)), ARF_RND_DOWN), m, z, prec);
        return;
    }

    /* todo: should be a test of whether n contains any negative integer */
    if (acb_contains_int(n) && !arb_is_nonnegative(acb_realref(n)))
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(t);
    acb_init(u);
    acb_init(v);

    acb_neg(t, n);
    acb_add_ui(u, m, 1, prec);
    acb_hypgeom_m(t, t, u, z, 1, prec);
    acb_add_ui(u, n, 1, prec);
    acb_rising(u, u, m, prec);
    acb_mul(res, t, u, prec);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

