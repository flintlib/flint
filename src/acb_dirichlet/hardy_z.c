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
acb_dirichlet_hardy_z(acb_ptr res, const acb_t t,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    slong len, slong prec)
{
    acb_ptr v, w;
    slong k;
    int is_real;

    if (len <= 0)
        return;

    /* use reflection formula -- todo for other characters */
    if ((G == NULL || G->q == 1) && arf_sgn(arb_midref(acb_imagref(t))) > 0)
    {
        acb_neg(res, t);
        acb_dirichlet_hardy_z(res, res, G, chi, len, prec);
        for (k = 1; k < len; k += 2)
            acb_neg(res + k, res + k);
        return;
    }

    v = _acb_vec_init(len);
    w = _acb_vec_init(len);

    is_real = acb_is_real(t);

    /* v = exp(i theta(t+x)) */
    acb_dirichlet_hardy_theta(v, t, G, chi, len, prec);
    _acb_vec_scalar_mul_onei(v, v, len);
    _acb_poly_exp_series(v, v, len, len, prec);

    /* w = L(1/2 + i (t+x)) */
    acb_one(w);
    acb_mul_2exp_si(w, w, -1);
    arb_sub(acb_realref(w), acb_realref(w), acb_imagref(t), prec);
    arb_set(acb_imagref(w), acb_realref(t));
    acb_dirichlet_l_jet(w, w, G, chi, 0, len, prec);
    for (k = 0; k < len; k++)
    {
        if (k % 4 == 1)
            acb_mul_onei(w + k, w + k);
        else if (k % 4 == 2)
            acb_neg(w + k, w + k);
        else if (k % 4 == 3)
            acb_div_onei(w + k, w + k);
    }

    _acb_poly_mullow(res, v, len, w, len, len, prec);

    if (is_real)
        for (k = 0; k < len; k++)
            arb_zero(acb_imagref(res + k));

    _acb_vec_clear(v, len);
    _acb_vec_clear(w, len);
}

