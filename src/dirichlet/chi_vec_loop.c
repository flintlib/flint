/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

static void
dirichlet_exponents_char(ulong * expo, const dirichlet_group_t G, const dirichlet_char_t chi, ulong order)
{
    slong k;
    ulong factor = G->expo / order;
    for (k = 0; k < G->num; k++)
        /* no overflow: log[k] < phi[k] and G->expo = phi[k] * PHI[k] */
        expo[k] = (chi->log[k] * G->PHI[k]) / factor;
}

/* loop over whole group */
void
dirichlet_chi_vec_loop_order(ulong * v, const dirichlet_group_t G, const dirichlet_char_t chi, ulong order, slong nv)
{
    int j;
    ulong t;
    slong k;
    ulong expo[MAX_FACTORS];
    dirichlet_char_t x;
    nmod_t o;

    dirichlet_char_init(x, G);
    dirichlet_char_one(x, G);

    dirichlet_exponents_char(expo, G, chi, order);
    nmod_init(&o, order);

    for (k = 0; k < nv; k++)
        v[k] = DIRICHLET_CHI_NULL;

    t = v[1] = 0;

    while ( (j = dirichlet_char_next(x, G)) >= 0 )
    {
        /* exponents were modified up to j */
        for (k = G->num - 1; k >= j; k--)
            t = nmod_add(t, expo[k], o);

        if (x->n < nv)
            v[x->n] = t;
    }

    /* fix result outside primes */
    /* dirichlet_vec_set_null(v, G, nv);*/
    /* copy outside modulus */

    for (k = G->q; k < nv ; k++ )
        v[k] = v[k - G->q];

    dirichlet_char_clear(x);
}

void
dirichlet_chi_vec_loop(ulong * v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv)
{
    dirichlet_chi_vec_loop_order(v, G, chi, G->expo, nv);
}
