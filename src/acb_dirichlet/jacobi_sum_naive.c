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
acb_dirichlet_jacobi_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)
{

    ulong k1, k2, m1, m2, g, e, m;
    ulong * v1, * v2;
    slong *v;
    nmod_t expo;
    acb_t z;

    v1 = flint_malloc(G->q * sizeof(ulong));
    v2 = flint_malloc(G->q * sizeof(ulong));

    dirichlet_vec_set_null(v1, G, G->q);
    dirichlet_chi_vec_loop(v1, G, chi1, G->q);

    dirichlet_vec_set_null(v2, G, G->q);
    dirichlet_chi_vec_loop(v2, G, chi2, G->q);

    nmod_init(&expo, G->expo);
    m1 = dirichlet_order_char(G, chi1);
    m2 = dirichlet_order_char(G, chi2);
    g = m1 * m2 / n_gcd(m1, m2);
    m = G->expo / g;

    v = flint_malloc(g * sizeof(slong));

    for (e = 0; e < g; e++)
        v[e] = 0;

    for (k1 = 2, k2 = G->q - 1; k2 > 1; k1++, k2--)
    {
        if (v1[k1] == DIRICHLET_CHI_NULL ||
            v2[k2] == DIRICHLET_CHI_NULL)
            continue;
        e = nmod_add(v1[k1], v2[k2], expo) / m;
        v[e]++;
    }

    acb_init(z);
    acb_unit_root(z, g, prec);
    acb_dirichlet_si_poly_evaluate(res, v, g, z, prec);

    acb_clear(z);
    flint_free(v);
    flint_free(v2);
    flint_free(v1);
}
