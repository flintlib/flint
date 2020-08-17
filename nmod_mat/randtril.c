/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"

void
nmod_mat_randtril(nmod_mat_t mat, flint_rand_t state, int unit)
{
    slong i, j;

    for (i = 0; i < mat->r; i++)
    {
        for (j = 0; j < mat->c; j++)
        {
            if (j < i)
            {
                nmod_mat_entry(mat, i, j) = n_randlimb(state) % (mat->mod.n);
            }
            else if (i == j)
            {
                nmod_mat_entry(mat, i, j) = n_randlimb(state) % (mat->mod.n);
                if (unit || nmod_mat_entry(mat, i, j) == UWORD(0))
                    nmod_mat_entry(mat, i, j) = UWORD(1);
            }
            else
            {
                nmod_mat_entry(mat, i, j) = UWORD(0);
            }
        }
    }
}
