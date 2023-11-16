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
acb_dirichlet_gauss_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    acb_t z;
    acb_ptr v;

    v = _acb_vec_init(G->q);
    acb_dirichlet_chi_vec(v, G, chi, G->q, prec);

    acb_init(z);
    acb_unit_root(z, G->q, prec);

    _acb_poly_evaluate(res, v, G->q, z, prec);

    acb_clear(z);
    _acb_vec_clear(v, G->q);
}
