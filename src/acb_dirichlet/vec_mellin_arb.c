/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_vec_mellin_arb(acb_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, const arb_t t, slong n, slong prec)
{
    slong k;
    arb_t tk, xt, stk, st;
    acb_ptr a;
    mag_t e;
    a = _acb_vec_init(len);
    acb_dirichlet_chi_vec(a, G, chi, len, prec);
    if (dirichlet_parity_char(G, chi))
    {
        for (k = 2; k < len; k++)
            acb_mul_si(a + k, a + k, k, prec);
    }
    arb_init(tk);
    arb_init(xt);
    arb_init(st);
    arb_init(stk);
    mag_init(e);

    arb_sqrt(st, t, prec);
    arb_one(tk);
    arb_one(stk);
    for (k = 0; k < n; k++)
    {
        _acb_dirichlet_theta_argument_at_arb(xt, G->q, tk, prec);
        mag_tail_kexpk2_arb(e, xt, len);
        arb_neg(xt, xt);
        arb_exp(xt, xt, prec);
        /* TODO: reduce len */
        acb_dirichlet_qseries_arb(res + k, a, xt, len, prec);
        acb_add_error_mag(res + k, e);
        acb_mul_arb(res + k, res + k, stk, prec);
        arb_mul(tk, tk, t, prec);
        arb_mul(stk, stk, st, prec);
    }
    mag_clear(e);
    arb_clear(xt);
    arb_clear(tk);
    arb_clear(stk);
    arb_clear(st);
    _acb_vec_clear(a, len);
}
