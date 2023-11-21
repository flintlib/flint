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
dirichlet_chi_vec_evenpart(ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, ulong order, slong nv)
{
    ulong mult = G->expo / order;
    if (G->neven >= 1 && chi->log[0])
    {
        ulong x, c3 = G->PHI[0] / mult;
        for (x = 3; x < nv; x += 4)
            v[x] = c3;
    }
    if (G->neven == 2 && chi->log[1])
    {
        ulong x, g, vx, xp, c4;
        nmod_t pe, o;
        nmod_init(&o, order);

        pe = G->P[1].pe;
        g = G->P[1].g;

        vx = c4 = (chi->log[1] * G->PHI[1]) / mult;
        for (x = g; x > 1;)
        {

            for (xp = x; xp < nv; xp += pe.n)
                v[xp] = nmod_add(v[xp], vx, o);

            for (xp = pe.n - x; xp < nv; xp += pe.n)
                v[xp] = nmod_add(v[xp], vx, o);

            x = nmod_mul(x, g, pe);
            vx = nmod_add(vx,  c4, o);
        }
    }
}

/* loop over primary components */
void
dirichlet_chi_vec_primeloop_order(ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, ulong order, slong nv)
{
    slong k, l;
    ulong mult = G->expo / order;
    nmod_t o;
    nmod_init(&o, order);

    for (k = 0; k < nv; k++)
        v[k] = 0;

    if (G->neven)
        dirichlet_chi_vec_evenpart(v, G, chi, order, nv);

    for (l = G->neven; l < G->num; l++)
    {
        dirichlet_prime_group_struct P = G->P[l];

        /* FIXME: there may be some precomputed dlog in P if needed */
        if (P.dlog == NULL)
            dlog_vec_add(v, nv, P.g, (chi->log[l] * G->PHI[l]) / mult, P.pe, P.phi.n, o);
        else
            dlog_vec_add_precomp(v, nv, P.dlog, P.g, (chi->log[l] * G->PHI[l]) / mult, P.pe, P.phi.n, o);

    }
    dirichlet_vec_set_null(v, G, nv);
}

void
dirichlet_chi_vec_primeloop(ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv)
{
    dirichlet_chi_vec_primeloop_order(v, G, chi, G->expo, nv);
}
