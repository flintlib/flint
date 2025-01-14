/*
    Copyright (C) 2010, 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

void
nmod_mat_randfull(nmod_mat_t mat, flint_rand_t state)
{
    slong i, j;

    for (i = 0; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) = FLINT_MAX(1, n_randint(state, mat->mod.n));
}

void
nmod_mat_randops(nmod_mat_t mat, flint_rand_t state, slong count)
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
                    nmod_mat_entry(mat, j, k) = nmod_add(nmod_mat_entry(mat, j, k),
                        nmod_mat_entry(mat, i, k), mat->mod);
            else
                for (k = 0; k < n; k++)
                    nmod_mat_entry(mat, j, k) = nmod_sub(nmod_mat_entry(mat, j, k),
                        nmod_mat_entry(mat, i, k), mat->mod);
        }
        else
        {
            if ((i = n_randint(state, n)) == (j = n_randint(state, n)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < m; k++)
                    nmod_mat_entry(mat, k, j) = nmod_add(nmod_mat_entry(mat, k, j),
                        nmod_mat_entry(mat, k, i), mat->mod);

            else
                for (k = 0; k < m; k++)
                    nmod_mat_entry(mat, k, j) = nmod_sub(nmod_mat_entry(mat, k, j),
                        nmod_mat_entry(mat, k, i), mat->mod);
        }
    }
}

int
nmod_mat_randpermdiag(nmod_mat_t mat, flint_rand_t state,
                            nn_srcptr diag, slong n)
{
    int parity;
    slong i;
    slong * rows;
    slong * cols;

    rows = _perm_init(mat->r);
    cols = _perm_init(mat->c);

    parity = _perm_randtest(rows, mat->r, state);
    parity ^= _perm_randtest(cols, mat->c, state);

    nmod_mat_zero(mat);
    for (i = 0; i < n; i++)
        nmod_mat_entry(mat, rows[i], cols[i]) = diag[i];

    _perm_clear(rows);
    _perm_clear(cols);

    return parity;
}

void
nmod_mat_randrank(nmod_mat_t mat, flint_rand_t state, slong rank)
{
    slong i;
    ulong * diag;

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

void
nmod_mat_randtest(nmod_mat_t mat, flint_rand_t state)
{
    _nmod_vec_randtest(mat->entries, state, mat->r * mat->c, mat->mod);
}

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

void
nmod_mat_randtriu(nmod_mat_t mat, flint_rand_t state, int unit)
{
    slong i, j;

    for (i = 0; i < mat->r; i++)
    {
        for (j = 0; j < mat->c; j++)
        {
            if (j > i)
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
