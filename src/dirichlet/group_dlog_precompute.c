/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

void
dirichlet_prime_group_dlog_precompute(dirichlet_prime_group_struct * P, ulong num)
{
    P->dlog = flint_malloc(sizeof(dlog_precomp_t));
    /* even if e = 1 use modpe struct, more flexible if reused in bigger group */
    dlog_precomp_modpe_init(P->dlog, P->g, P->p, P->e, P->pe.n, num);
}

void
dirichlet_group_dlog_precompute(dirichlet_group_t G, ulong num)
{
    slong k;
    for (k = 0; k < G->num; k++)
    {
        if (G->P[k].dlog == NULL)
            dirichlet_prime_group_dlog_precompute(&G->P[k], num);
    }
}

void
dirichlet_group_dlog_clear(dirichlet_group_t G)
{
    slong k;
    for (k = 0; k < G->num; k++)
    {
        if (G->P[k].dlog != NULL)
        {
            dlog_precomp_clear(G->P[k].dlog);
            flint_free(G->P[k].dlog);
            G->P[k].dlog = NULL;
        }
    }
}
