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
use_recurrence(const acb_t n, const acb_t a, const acb_t b, slong prec)
{
    if (!acb_is_int(n) || !arb_is_nonnegative(acb_realref(n)))
        return 0;

    if (arf_cmpabs_ui(arb_midref(acb_realref(n)), prec) > 0)
        return 0;

    if (arb_is_nonnegative(acb_realref(a)) ||
        arf_get_d(arb_midref(acb_realref(a)), ARF_RND_DOWN) > -0.9)
        return 0;

    return 1;
}

void
acb_hypgeom_jacobi_p_ui_direct(acb_t res, ulong n,
    const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    acb_ptr terms;
    acb_t t, u, v;
    slong k;

    terms = _acb_vec_init(n + 1);
    acb_init(t);
    acb_init(u);
    acb_init(v);

    acb_one(terms);
    acb_add_ui(u, z, 1, prec);

    for (k = 1; k <= n; k++)
    {
        acb_add_ui(t, a, n + 1 - k, prec);
        acb_mul(t, t, u, prec);
        acb_div_ui(t, t, 2 * k, prec);
        acb_mul(terms + k, terms + k - 1, t, prec);
    }

    acb_sub_ui(u, z, 1, prec);
    acb_one(v);

    for (k = 1; k <= n; k++)
    {
        acb_add_ui(t, b, n + 1 - k, prec);
        acb_mul(t, t, u, prec);
        acb_div_ui(t, t, 2 * k, prec);
        acb_mul(v, v, t, prec);
        acb_mul(terms + n - k, terms + n - k, v, prec);
    }

    acb_set(res, terms);
    for (k = 1; k <= n; k++)
        acb_add(res, res, terms + k, prec);

    _acb_vec_clear(terms, n + 1);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
acb_hypgeom_jacobi_p(acb_t res, const acb_t n, const acb_t a, const acb_t b,
    const acb_t z, slong prec)
{
    acb_t t, u, v, w;

    if (use_recurrence(n, a, b, prec))
    {
        acb_hypgeom_jacobi_p_ui_direct(res,
            arf_get_si(arb_midref(acb_realref(n)), ARF_RND_DOWN), a, b, z, prec);
        return;
    }

    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_init(w);

    acb_neg(t, n);
    acb_add_ui(v, a, 1, prec);
    acb_add(u, n, v, prec);
    acb_add(u, u, b, prec);
    acb_sub_ui(w, z, 1, prec);
    acb_mul_2exp_si(w, w, -1);
    acb_neg(w, w);

    acb_hypgeom_2f1(w, t, u, v, w, 0, prec);

    acb_rising(t, v, n, prec);
    acb_mul(w, w, t, prec);

    acb_add_ui(t, n, 1, prec);
    acb_rgamma(t, t, prec);
    acb_mul(w, w, t, prec);

    acb_set(res, w);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    acb_clear(w);
}

