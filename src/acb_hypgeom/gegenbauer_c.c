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

    if (arb_is_nonnegative(acb_realref(m)))
        return 0;

    return 1;
}

void
acb_hypgeom_gegenbauer_c_ui_recurrence(acb_t res, ulong n, const acb_t m,
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
        acb_mul(res, m, z, prec);
        acb_mul_2exp_si(res, res, 1);
        return;
    }

    acb_init(t);
    acb_init(u);
    acb_init(v);

    acb_one(t);
    acb_mul(u, m, z, prec);
    acb_mul_2exp_si(u, u, 1);

    for (k = 2; k <= n; k++)
    {
        acb_mul_2exp_si(v, m, 1);
        acb_add_ui(v, v, k - 2, prec);
        acb_mul(t, t, v, prec);

        acb_add_ui(v, m, k - 1, prec);
        acb_mul(v, v, z, prec);
        acb_mul_2exp_si(v, v, 1);
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
acb_hypgeom_gegenbauer_c(acb_t res, const acb_t n, const acb_t m,
    const acb_t z, slong prec)
{
    acb_t a, b, c, t;

    if (use_recurrence(n, m, prec))
    {
        acb_hypgeom_gegenbauer_c_ui_recurrence(res,
            arf_get_si(arb_midref(acb_realref(n)), ARF_RND_DOWN), m, z, prec);
        return;
    }

    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(t);

    acb_neg(a, n);
    acb_mul_2exp_si(b, m, 1);
    acb_add(b, b, n, prec);
    acb_one(c);
    acb_mul_2exp_si(c, c, -1);
    acb_add(c, c, m, prec);
    acb_sub_ui(t, z, 1, prec);
    acb_mul_2exp_si(t, t, -1);
    acb_neg(t, t);

    acb_hypgeom_2f1(t, a, b, c, t, 0, prec);

    acb_mul_2exp_si(b, m, 1);
    acb_rising(b, b, n, prec);
    acb_mul(t, t, b, prec);

    acb_add_ui(b, n, 1, prec);
    acb_rgamma(b, b, prec);
    acb_mul(t, t, b, prec);

    acb_set(res, t);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(t);
}

