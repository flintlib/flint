/*
    Copyright (C) 2010, 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2017 Luca De Feo
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_mod_mat.h"

void fmpz_mod_mat_randrank(fmpz_mod_mat_t mat, flint_rand_t state, slong rank)
{
    slong i;
    fmpz * diag;

    if (rank < 0 || rank > fmpz_mod_mat_nrows(mat) || rank > fmpz_mod_mat_ncols(mat))
        flint_throw(FLINT_ERROR, "Impossible rank in %s\n", __func__);

    diag = _fmpz_vec_init(rank);

    for (i = 0; i < rank; i++)
    {
        fmpz_randm(diag + i, state, mat->mod);
        if (fmpz_is_zero(diag + i))
            fmpz_one(diag + i);
    }

    fmpz_mat_randpermdiag(mat->mat, state, diag, rank);

    _fmpz_vec_clear(diag, rank);
}

void fmpz_mod_mat_randtest(fmpz_mod_mat_t mat, flint_rand_t state)
{
    fmpz_mat_randtest(mat->mat, state, fmpz_bits(mat->mod));
    _fmpz_mod_mat_reduce(mat);
}

void fmpz_mod_mat_randtril(fmpz_mod_mat_t mat, flint_rand_t state, int unit)
{
    fmpz* e;
    slong i, j;

    for (i = 0; i < fmpz_mod_mat_nrows(mat); i++)
    {
        for (j = 0; j < fmpz_mod_mat_ncols(mat); j++)
        {
            e = fmpz_mod_mat_entry(mat, i, j);
            if (j < i)
            {
                fmpz_randm(e, state, mat->mod);
            }
            else if (i == j)
            {
                fmpz_randm(e, state, mat->mod);
                if (unit || fmpz_is_zero(e))
                    fmpz_one(e);
            }
            else
            {
                fmpz_zero(e);
            }
        }
    }
}

void fmpz_mod_mat_randtriu(fmpz_mod_mat_t mat, flint_rand_t state, int unit)
{
    fmpz* e;
    slong i, j;

    for (i = 0; i < fmpz_mod_mat_nrows(mat); i++)
    {
        for (j = 0; j < fmpz_mod_mat_ncols(mat); j++)
        {
            e = fmpz_mod_mat_entry(mat, i, j);
            if (j > i)
            {
                fmpz_randm(e, state, mat->mod);
            }
            else if (i == j)
            {
                fmpz_randm(e, state, mat->mod);
                if (unit || fmpz_is_zero(e))
                    fmpz_one(e);
            }
            else
            {
                fmpz_zero(e);
            }
        }
    }
}

void fmpz_mod_mat_randops(fmpz_mod_mat_t mat, slong count, flint_rand_t state)
{
    fmpz_mat_randops(mat->mat, state, count);
    _fmpz_mod_mat_reduce(mat);
}
