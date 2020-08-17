/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod_vec.h"



void
nmod_mat_randops(nmod_mat_t mat, slong count, flint_rand_t state)
{
    slong c, i, j, k;
    slong m = mat->r;
    slong n = mat->c;

    if (mat->r == 0 || mat->c == 0)
        return;

    for (c = 0; c < count; c++)
    {
        if (n_randint(state, 2))
        {
            if ((i = n_randint(state, m)) == (j = n_randint(state, m)))
                continue;

            if (n_randint(state, 2))
                for (k = 0; k < n; k++)
                    mat->rows[j][k] = nmod_add(mat->rows[j][k],
                        mat->rows[i][k], mat->mod);
            else
                for (k = 0; k < n; k++)
                    mat->rows[j][k] = nmod_sub(mat->rows[j][k],
                        mat->rows[i][k], mat->mod);
        }
        else
        {
            if ((i = n_randint(state, n)) == (j = n_randint(state, n)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < m; k++)
                    mat->rows[k][j] = nmod_add(mat->rows[k][j],
                        mat->rows[k][i], mat->mod);
            else
                for (k = 0; k < m; k++)
                    mat->rows[k][j] = nmod_sub(mat->rows[k][j],
                        mat->rows[k][i], mat->mod);
        }
    }
}
