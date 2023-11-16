/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"
#include "acb_dirichlet.h"

/* dft, lexicographic conrey indexing, array size G->phi_q */
void
acb_dirichlet_dft_index(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec)
{
    if (G->phi_q == 1)
    {
        acb_set(w, v);
    }
    else
    {
        slong k, l, * cyc;


        cyc = flint_malloc(G->num * sizeof(slong));
        for (k = 0, l = G->num - 1; l >= 0; k++, l--)
            cyc[k] = G->P[k].phi.n;

        acb_dft_prod(w, v, cyc, G->num, prec);
        flint_free(cyc);
    }
}

/* dft, number indexing, array size G->q */
void
acb_dirichlet_dft(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec)
{
    ulong i, len;
    acb_ptr t1, t2;
    dirichlet_char_t x;

    len = G->phi_q;
    t1 = flint_malloc(len * sizeof(acb_struct));

    dirichlet_char_init(x, G);
    dirichlet_char_one(x, G);
    for (i = 0; i < len; i++)
    {
        t1[i] = v[x->n];
        dirichlet_char_next(x, G);
    };

    t2 = _acb_vec_init(len);
    acb_dirichlet_dft_index(t2, t1, G, prec);

    dirichlet_char_one(x, G);
    for (i = 0; i < len; i++)
    {
        acb_set(w + x->n, t2 + i);
        dirichlet_char_next(x, G);
    };

    _acb_vec_clear(t2, len);
    dirichlet_char_clear(x);
    flint_free(t1);
}
