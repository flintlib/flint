/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "nmod_vec.h"

void
nmod_mat_randrank(nmod_mat_t mat, flint_rand_t state, slong rank)
{
    slong i;
    mp_limb_t * diag;

    if (rank < 0 || rank > mat->r || rank > mat->c)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_mat_randrank). Impossible rank.\n");
    }

    diag = _nmod_vec_init(rank);

    if (mat->mod.n != 1)
    {
        for (i = 0; i < rank; i++)
            diag[i] = 1 + n_randint(state, mat->mod.n - 1);
    } else
    {
        for (i = 0; i < rank; i++)
            diag[i] = 0;
    }

    nmod_mat_randpermdiag(mat, state, diag, rank);

    _nmod_vec_clear(diag);
}
