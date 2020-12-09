/*
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

void
nmod_mat_scalar_addmul_ui(nmod_mat_t dest, const nmod_mat_t X,
                                          const nmod_mat_t Y, const mp_limb_t b)
{
    slong i, j;

    if (b == UWORD(0))
    {
        if (dest != X)
            nmod_mat_set(dest, X);
        return;
    }

    for (i = 0; i < X->r; i++)
    {
        for (j = 0; j < X->c; j++)
        {
            nmod_mat_entry(dest, i, j) =
                n_addmod(nmod_mat_entry(X, i, j),
                    n_mulmod2_preinv(nmod_mat_entry(Y, i, j), b, Y->mod.n,
                                                        Y->mod.ninv), X->mod.n);
        }
    }
}
