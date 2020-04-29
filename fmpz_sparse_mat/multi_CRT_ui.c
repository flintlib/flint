/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"

void
fmpz_sparse_mat_multi_CRT_ui_precomp(fmpz_sparse_mat_t M,
    nmod_sparse_mat_struct * residues, slong nres,
    const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign)
{
    slong i, j;
    nmod_sparse_vec_struct *vres;

    vres = flint_malloc(nres * sizeof(*vres));

    for (i = 0; i < M->r; i++)
    {
        for (j = 0; j < nres; ++j) vres[j] = residues[j].rows[i];
        fmpz_sparse_vec_multi_CRT_ui_precomp(&M->rows[i], vres, nres, comb, temp, sign);
    }
    flint_free(vres);
}

void
fmpz_sparse_mat_multi_CRT_ui(fmpz_sparse_mat_t mat, nmod_sparse_mat_struct * residues, slong nres, int sign)
{
    slong i;
    mp_limb_t *primes;
    fmpz_comb_t comb;
    fmpz_comb_temp_t temp;

    primes = flint_malloc(nres * sizeof(*primes));
    for (i = 0; i < nres; i++) primes[i] = residues[i].mod.n;
    fmpz_comb_init(comb, primes, nres);
    fmpz_comb_temp_init(temp, comb);

    fmpz_sparse_mat_multi_CRT_ui_precomp(mat, residues, nres, comb, temp, sign);

    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(temp);
    flint_free(primes);
}
