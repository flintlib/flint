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
acb_dirichlet_root_number_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    arb_t x;
    acb_t eps;

    arb_init(x);
    arb_one(x);
    acb_dirichlet_theta_arb(res, G, chi, x, prec);
    acb_init(eps);
    acb_conj(eps, res);
    acb_div(res, res, eps, prec);

    if (dirichlet_char_is_real(G, chi))
        arb_zero(acb_imagref(res));

    arb_clear(x);
    acb_clear(eps);
}

void
acb_dirichlet_root_number(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    if (dirichlet_conductor_char(G, chi) < G->q)
    {
        flint_throw(FLINT_ERROR, "root number: need primitive character\n");
    }
    else if (G->num > 1)
    {
        acb_t iq;
        acb_init(iq);
        acb_dirichlet_gauss_sum_order2(iq, G, chi, prec);
        acb_dirichlet_gauss_sum(res, G, chi, prec);
        acb_div(res, res, iq, prec);
        acb_clear(iq);
    }
    else
    {
        acb_dirichlet_root_number_theta(res, G, chi, prec);
    }

    if (dirichlet_char_is_real(G, chi))
        arb_zero(acb_imagref(res));
}
