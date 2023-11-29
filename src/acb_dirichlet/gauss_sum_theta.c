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
acb_dirichlet_gauss_sum_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    ulong cond = dirichlet_conductor_char(G, chi);
    if (cond < G->q || (G->q == 300 && (chi->n == 71 || chi->n == 131))
            || (G->q == 600 && (chi->n == 11 || chi->n == 491)))
    {
        flint_throw(FLINT_ERROR, "(gauss_sum_theta): non available for non primitive character"
                "or exceptional characters chi_300(71,.), chi_300(131,.), "
                "chi_600(11,.) and chi_600(491,.)\n");
    }
    else
    {
        acb_t iq;
        acb_init(iq);
        acb_dirichlet_gauss_sum_order2(iq, G, chi, prec);
        acb_dirichlet_root_number_theta(res, G, chi, prec);
        acb_mul(res, res, iq, prec);
        acb_clear(iq);
    }
}
