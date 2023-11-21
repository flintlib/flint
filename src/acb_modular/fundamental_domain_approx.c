/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_modular.h"

static int
good_enough(const acb_t z, const arf_t one_minus_eps, slong prec)
{
    arf_t m;
    int res;

    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -1) > 0)
        return 0;

    if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) >= 0)
        return 1;

    arf_init(m);
    arf_mul(m, arb_midref(acb_realref(z)), arb_midref(acb_realref(z)), prec, ARF_RND_DOWN);
    arf_addmul(m, arb_midref(acb_imagref(z)), arb_midref(acb_imagref(z)), prec, ARF_RND_DOWN);
    res = (arf_cmp(m, one_minus_eps) >= 0);
    arf_clear(m);

    return res;
}

void
acb_modular_fundamental_domain_approx(acb_t w, psl2z_t g, const acb_t z,
        const arf_t one_minus_eps, slong prec)
{
    acb_t t;

    psl2z_one(g);

    /* we must be in the upper half-plane */
    if (!acb_is_finite(z) || arf_sgn(arb_midref(acb_imagref(z))) <= 0)
    {
        acb_set(w, z);
        return;
    }

    /* too large real-value shift */
    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), prec) > 0)
    {
        acb_set(w, z);
        return;
    }

    /* y >= 1: just shift x */
    if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) >= 0)
    {
        arf_get_fmpz(&g->b, arb_midref(acb_realref(z)), ARF_RND_NEAR);
        acb_sub_fmpz(w, z, &g->b, prec);
        fmpz_neg(&g->b, &g->b);
        return;
    }

    acb_init(t);

    /* try using doubles */
    if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -40) > 0 &&
        arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 40) < 0)
    {
        double zred, zimd;

        zred = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
        zimd = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

        acb_modular_fundamental_domain_approx_d(g, zred, zimd,
            arf_get_d(one_minus_eps, ARF_RND_UP));
        acb_modular_transform(t, g, z, prec);

        if (good_enough(t, one_minus_eps, prec))
        {
            acb_swap(w, t);
            acb_clear(t);
            return;
        }
    }

    /* try with full precision */
    acb_modular_fundamental_domain_approx_arf(g,
        arb_midref(acb_realref(z)), arb_midref(acb_imagref(z)),
        one_minus_eps, prec);

    acb_modular_transform(t, g, z, prec);
    acb_swap(w, t);

    acb_clear(t);
}

