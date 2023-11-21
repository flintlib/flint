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
acb_dirichlet_chi_vec(acb_ptr v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv, slong prec)
{
    slong k;
    ulong * a, order;
    acb_dirichlet_roots_t t;

    a = flint_malloc(nv * sizeof(ulong));
    order = dirichlet_order_char(G, chi);
    dirichlet_chi_vec_order(a, G, chi, order, nv);

    acb_dirichlet_roots_init(t, order, nv, prec);

    acb_zero(v + 0);
    for (k = 0; k < nv; k++)
    {
        if (a[k] != DIRICHLET_CHI_NULL)
            acb_dirichlet_root(v + k, t, a[k], prec);
        else
            acb_zero(v + k);
    }

    acb_dirichlet_roots_clear(t);
    flint_free(a);
}
