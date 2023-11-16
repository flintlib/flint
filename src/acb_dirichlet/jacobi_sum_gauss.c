/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* should use only for prime power modulus */
void
acb_dirichlet_jacobi_sum_gauss(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)
{
    /* J_q(a,b)G_q(ab) = G_q(a)G_q(b) */
    acb_t tmp;
    dirichlet_char_t chi12;

    dirichlet_char_init(chi12, G);
    dirichlet_char_mul(chi12, G, chi1, chi2);

    acb_init(tmp);

    acb_dirichlet_gauss_sum(res, G, chi1, prec);
    if (chi2->n == chi1->n)
        acb_set(tmp, res);
    else
        acb_dirichlet_gauss_sum(tmp, G, chi2, prec);
    acb_mul(res, res, tmp, prec);
    acb_dirichlet_gauss_sum(tmp, G, chi12, prec);
    acb_div(res, res, tmp, prec);

    dirichlet_char_clear(chi12);
    acb_clear(tmp);
}
