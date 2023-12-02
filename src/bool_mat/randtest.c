/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "bool_mat.h"

void
bool_mat_randtest(bool_mat_t mat, flint_rand_t state)
{
    slong i, j;
    mp_limb_t density;

    density = n_randint(state, 101);
    for (i = 0; i < bool_mat_nrows(mat); i++)
        for (j = 0; j < bool_mat_ncols(mat); j++)
            bool_mat_set_entry(mat, i, j, n_randint(state, 100) < density);
}

void
bool_mat_randtest_diagonal(bool_mat_t mat, flint_rand_t state)
{
    slong n, i;
    slong density;

    n = FLINT_MIN(bool_mat_nrows(mat), bool_mat_ncols(mat));

    density = n_randint(state, 101);
    bool_mat_zero(mat);
    for (i = 0; i < n; i++)
        bool_mat_set_entry(mat, i, i, n_randint(state, 100) < density);
}

void
bool_mat_randtest_nilpotent(bool_mat_t mat, flint_rand_t state)
{
    slong n, i, j;
    slong density;

    if (!bool_mat_is_square(mat))
    {
        flint_throw(FLINT_ERROR, "bool_mat_randtest_nilpotent: "
                     "a square matrix is required!\n");
    }

    if (bool_mat_is_empty(mat))
    {
        flint_throw(FLINT_ERROR, "bool_mat_randtest_nilpotent: "
                     "a non-empty matrix is required!\n");
    }

    n = bool_mat_nrows(mat);

    if (n == 1)
    {
        bool_mat_set_entry(mat, 0, 0, 0);
        return;
    }

    /* sample a strictly triangular matrix */
    density = n_randint(state, 101);
    bool_mat_zero(mat);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            bool_mat_set_entry(mat, i, j, n_randint(state, 100) < density);
        }
    }

    /* permute rows and columns */
    {
        slong *p;
        bool_mat_t A;

        bool_mat_init(A, n, n);
        bool_mat_set(A, mat);

        p = flint_malloc(n * sizeof(slong));
        _perm_randtest(p, n, state);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                bool_mat_set_entry(
                        mat, p[i], p[j], bool_mat_get_entry(A, i, j));
            }
        }

        flint_free(p);
        bool_mat_clear(A);
    }
}
