/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "acb_poly.h"

void
_acb_dirichlet_theta_arb_smallorder(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const arb_t xt, slong len, slong prec)
{
    ulong order;
    ulong * a;
    int parity;
    acb_dirichlet_roots_t z;

    parity = dirichlet_parity_char(G, chi);
    order = dirichlet_order_char(G, chi);
    a = flint_malloc(len * sizeof(ulong));
    dirichlet_chi_vec_order(a, G, chi, order, len);

    acb_dirichlet_roots_init(z, order, 0, prec);
    acb_dirichlet_qseries_arb_powers_smallorder(res, xt, parity, a, z, len, prec);
    acb_dirichlet_roots_clear(z);

    flint_free(a);
}

void
_acb_dirichlet_theta_arb_series(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const arb_t xt, slong len, slong prec)
{
    acb_ptr a;
    a = _acb_vec_init(len);
    acb_dirichlet_chi_vec(a, G, chi, len, prec);
    if (dirichlet_parity_char(G, chi))
    {
        slong k;
        for (k = 2; k < len; k++)
            acb_mul_si(a + k, a + k, k, prec);
    }
    acb_dirichlet_qseries_arb(res, a, xt, len, prec);
    _acb_vec_clear(a, len);
}

void
_acb_dirichlet_theta_arb_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const arb_t xt, slong len, slong prec)
{
    ulong order, * a;
    acb_dirichlet_roots_t z;
    int parity;

    parity = dirichlet_parity_char(G, chi);
    order = dirichlet_order_char(G, chi);
    a = flint_malloc(len * sizeof(ulong));
    dirichlet_chi_vec_order(a, G, chi, order, len);

    acb_dirichlet_roots_init(z, order, len, prec);

    acb_dirichlet_qseries_arb_powers_naive(res, xt, parity, a, z, len, prec);

    acb_dirichlet_roots_clear(z);
    flint_free(a);
}

void
acb_dirichlet_theta_arb(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const arb_t t, slong prec)
{
    slong len;
    ulong order;
    arb_t xt;
    mag_t e;

    len = acb_dirichlet_theta_length(G->q, t, prec);

    arb_init(xt);
    _acb_dirichlet_theta_argument_at_arb(xt, G->q, t, prec);

    mag_init(e);
    mag_tail_kexpk2_arb(e, xt, len);

    arb_neg(xt, xt);
    arb_exp(xt, xt, prec);

    /* TODO: tune this limit */
    order = dirichlet_order_char(G, chi);
    if (order < 30)
        _acb_dirichlet_theta_arb_smallorder(res, G, chi, xt, len, prec);
    else
        _acb_dirichlet_theta_arb_naive(res, G, chi, xt, len, prec);

    arb_add_error_mag(acb_realref(res), e);
    arb_add_error_mag(acb_imagref(res), e);

    mag_clear(e);

    acb_mul_2exp_si(res, res, 1);
    arb_clear(xt);
}
