/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

/* order of an element knowing the factorization of a multiple */
ulong
nmod_order_precomp(ulong a, nmod_t mod, ulong expo, n_factor_t fac)
{
    int k;
    ulong pe, ap, order = 1;
    for (k = 0; k < fac.num; k++)
    {
        pe = n_pow(fac.p[k], fac.exp[k]);
        ap = nmod_pow_ui(a, expo / pe, mod);
        while ( ap != 1)
        {
            ap = nmod_pow_ui(ap, fac.p[k], mod);
            order *= fac.p[k];
        }
    }
    return order;
}

ulong
dirichlet_order_ui(const dirichlet_group_t G, ulong a)
{
    n_factor_t fac;

    n_factor_init(&fac);
    n_factor(&fac, G->expo, 1);
    return nmod_order_precomp(a, G->mod, G->expo, fac);
}
