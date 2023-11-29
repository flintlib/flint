/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "acb_poly.h"

void
acb_dirichlet_hardy_theta(acb_ptr res, const acb_t t,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    slong len, slong prec)
{
    acb_struct y[2];
    arb_t c;
    slong k;
    int parity;
    ulong q;

    if (len <= 0)
        return;

    if (t == res)
    {
        acb_init(y);
        acb_set(y, t);
        acb_dirichlet_hardy_theta(res, y, G, chi, len, prec);
        acb_clear(y);
        return;
    }

    if (G == NULL)
    {
        parity = 0;
        q = 1;
    }
    else
    {
        parity = dirichlet_parity_char(G, chi);
        q = G->q;

        if (q != dirichlet_conductor_char(G, chi))
        {
            flint_throw(FLINT_ERROR, "hardy theta: need primitive character\n");
        }
    }

    if (!acb_is_finite(t))
    {
        _acb_vec_indeterminate(res, len);
        return;
    }

    acb_init(y + 0);
    acb_init(y + 1);
    arb_init(c);

    /* res = log gamma((s+parity)/2), s = 0.5+i(t+x) */
    acb_mul_onei(y, t);
    arb_set_d(c, 0.5 + parity);
    arb_add(acb_realref(y), acb_realref(y), c, prec);
    acb_mul_2exp_si(y, y, -1);
    acb_onei(y + 1);
    acb_mul_2exp_si(y + 1, y + 1, -1);
    _acb_poly_lgamma_series(res, y, FLINT_MIN(len, 2), len, prec);

    if (acb_is_real(t))
    {
        for (k = 0; k < len; k++)
        {
            arb_set(acb_realref(res + k), acb_imagref(res + k));
            arb_zero(acb_imagref(res + k));
        }
    }
    else
    {
        acb_ptr v = _acb_vec_init(len);

        /* v = log gamma(((1-s)+parity)/2), s = 0.5+i(t+x) */
        acb_div_onei(y, t);
        arb_set_d(c, 0.5 + parity);
        arb_add(acb_realref(y), acb_realref(y), c, prec);
        acb_mul_2exp_si(y, y, -1);
        acb_neg(y + 1, y + 1);
        _acb_poly_lgamma_series(v, y, FLINT_MIN(len, 2), len, prec);

        _acb_vec_sub(res, res, v, len, prec);
        for (k = 0; k < len; k++)
        {
            acb_div_onei(res + k, res + k);
            acb_mul_2exp_si(res + k, res + k, -1);
        }

        _acb_vec_clear(v, len);
    }

    /* (t+x) [-(1/2) log(pi/q)] */
    arb_const_pi(c, prec);
    arb_div_ui(c, c, q, prec);
    arb_log(c, c, prec);
    arb_mul_2exp_si(c, c, -1);
    acb_submul_arb(res, t, c, prec);
    if (len > 1)
        acb_sub_arb(res + 1, res + 1, c, prec);

    /* i log(eps) / 2 */
    if (q > 1)
    {
        acb_dirichlet_root_number(y, G, chi, prec);
        acb_arg(c, y, prec);
        arb_mul_2exp_si(c, c, -1);
        arb_sub(acb_realref(res), acb_realref(res), c, prec);
    }

    acb_clear(y + 0);
    acb_clear(y + 1);
    arb_clear(c);
}

