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

void
acb_hypgeom_hermite_h_ui_recurrence(acb_t res, ulong n, const acb_t z, slong prec)
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
        acb_set_round(res, z, prec);
        acb_mul_2exp_si(res, res, 1);
        return;
    }

    acb_init(t);
    acb_init(u);
    acb_init(v);

    acb_one(t);
    acb_mul_2exp_si(u, z, 1);

    for (k = 2; k <= n; k++)
    {
        acb_mul(v, u, z, prec);
        acb_submul_ui(v, t, k - 1, prec);
        acb_mul_2exp_si(v, v, 1);
        acb_swap(t, u);
        acb_swap(u, v);
    }

    acb_set(res, u);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
acb_hypgeom_hermite_h(acb_t res, const acb_t n, const acb_t z, slong prec)
{
    acb_t a, b, c, t, u, v;
    int use_asymp;

    if (acb_is_int(n) && arb_is_nonnegative(acb_realref(n)) &&
        (arf_cmpabs_ui(arb_midref(acb_realref(n)), 30) < 0))
    {
        acb_hypgeom_hermite_h_ui_recurrence(res,
            arf_get_si(arb_midref(acb_realref(n)), ARF_RND_DOWN), z, prec);
        return;
    }

    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(t);
    acb_init(u);
    acb_init(v);

    acb_mul(t, z, z, prec);

    use_asymp = arb_is_positive(acb_realref(z)) &&
        acb_hypgeom_u_use_asymp(t, prec);

    if (use_asymp)
    {
        acb_mul_2exp_si(a, n, -1);
        acb_neg(a, a);
        acb_one(b);
        acb_mul_2exp_si(b, b, -1);
        acb_hypgeom_u_asymp(u, a, b, t, -1, prec);
        acb_mul_2exp_si(t, z, 1);
        acb_pow(t, t, n, prec);
        acb_mul(u, u, t, prec);
        acb_set(res, u);
    }
    else
    {
        /* a = (1-n)/2 */
        acb_sub_ui(a, n, 1, prec);
        acb_neg(a, a);
        acb_mul_2exp_si(a, a, -1);
        /* c = -n/2 */
        acb_mul_2exp_si(c, n, -1);
        acb_neg(c, c);

        acb_rgamma(u, a, prec);

        if (!acb_is_zero(u))
        {
            acb_one(b);
            acb_mul_2exp_si(b, b, -1);
            acb_hypgeom_m(v, c, b, t, 0, prec);
            acb_mul(u, u, v, prec);
        }

        acb_rgamma(v, c, prec);

        if (!acb_is_zero(v))
        {
            acb_set_ui(b, 3);
            acb_mul_2exp_si(b, b, -1);
            acb_hypgeom_m(t, a, b, t, 0, prec);
            acb_mul_2exp_si(t, t, 1);
            acb_mul(t, t, z, prec);
            acb_submul(u, v, t, prec);
        }

        acb_set_ui(t, 2);
        acb_pow(t, t, n, prec);
        acb_mul(u, u, t, prec);
        arb_const_sqrt_pi(acb_realref(t), prec);
        acb_mul_arb(u, u, acb_realref(t), prec);

        acb_set(res, u);
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

