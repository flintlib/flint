/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_naive_term(acb_t res, acb_srcptr z, const acb_mat_t tau,
    slong* n, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_ptr x, y, v;
    arb_mat_t X, Y;
    acb_t dot;
    slong k;

    x = _arb_vec_init(g);
    y = _arb_vec_init(g);
    v = _arb_vec_init(g);
    arb_mat_init(X, g, g);
    arb_mat_init(Y, g, g);
    acb_init(dot);

    _acb_vec_get_real(x, z, g);
    _acb_vec_get_imag(y, z, g);
    acb_mat_get_real(X, tau);
    acb_mat_get_imag(Y, tau);    
    for (k = 0; k < g; k++)
    {
        arb_set_si(&v[k], n[k]);
    }

    acb_zero(res);
    arb_mat_bilinear_form(acb_realref(res), X, v, v, prec);
    arb_mat_bilinear_form(acb_imagref(res), Y, v, v, prec);
    arb_dot(acb_realref(dot), NULL, 0, v, 1, x, 1, g, prec);
    arb_dot(acb_imagref(dot), NULL, 0, v, 1, y, 1, g, prec);
    arb_mul_2exp_si(dot, dot, 1);
    acb_add(res, res, dot, prec);
    acb_exp_pi_i(res, res, prec);

    _arb_vec_clear(x, g);
    _arb_vec_clear(y, g);
    _arb_vec_clear(v, g);
    arb_mat_clear(X);
    arb_mat_clear(Y);
    acb_clear(dot);
}
