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
acb_hypgeom_chebyshev_u(acb_t res, const acb_t n, const acb_t z, slong prec)
{
    acb_t t, u;

    if (acb_is_int(n) &&
        arf_cmpabs_2exp_si(arb_midref(acb_realref(n)), FLINT_BITS - 1) < 0)
    {
        slong k = arf_get_si(arb_midref(acb_realref(n)), ARF_RND_DOWN);

        if (k >= 0)
        {
            acb_chebyshev_u_ui(res, k, z, prec);
        }
        else if (k == -1)
        {
            acb_zero(res);
        }
        else
        {
            acb_chebyshev_u_ui(res, -2-k, z, prec);
            acb_neg(res, res);
        }

        return;
    }

    if (acb_is_zero(z))
    {
        acb_mul_2exp_si(res, n, -1);
        acb_cos_pi(res, res, prec);
        return;
    }

    if (acb_is_one(z))
    {
        acb_add_ui(res, n, 1, prec);
        return;
    }

    acb_init(t);
    acb_init(u);

    acb_add_ui(u, n, 1, prec);
    acb_sub_ui(t, z, 1, prec);

    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(t)), -2 - prec / 10) < 0 &&
        arf_cmpabs_2exp_si(arb_midref(acb_imagref(t)), -2 - prec / 10) < 0)
    {
        acb_t a, b, c;

        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_neg(a, n);
        acb_add_ui(b, n, 2, prec);
        acb_set_ui(c, 3);
        acb_mul_2exp_si(c, c, -1);
        acb_neg(t, t);
        acb_mul_2exp_si(t, t, -1);
        acb_hypgeom_2f1(t, a, b, c, t, 0, prec);
        acb_mul(res, t, u, prec);

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }
    else
    {
        if (arb_is_positive(acb_realref(t)))
        {
            /* sinh((n+1) acosh(z)) / (sqrt(z-1) sqrt(z+1));
               can use one square root when strictly in the right half plane */
            acb_mul(t, z, z, prec);
            acb_sub_ui(t, t, 1, prec);
            acb_acosh(res, z, prec);
            acb_mul(res, res, u, prec);
            acb_sinh(res, res, prec);
            acb_rsqrt(t, t, prec);
            acb_mul(res, res, t, prec);
        }
        else
        {
            acb_mul(t, z, z, prec);
            acb_sub_ui(t, t, 1, prec);
            acb_acos(res, z, prec);
            acb_mul(res, res, u, prec);
            acb_sin(res, res, prec);
            acb_neg(t, t);
            acb_rsqrt(t, t, prec);
            acb_mul(res, res, t, prec);
        }
    }

    acb_clear(t);
    acb_clear(u);
}

