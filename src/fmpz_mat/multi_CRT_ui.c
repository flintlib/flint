/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_multi_CRT_ui_precomp(fmpz_mat_t mat,
    nmod_mat_t * const residues, slong nres,
    const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign)
{
    slong i, j, k;
    mp_ptr r;

    r = _nmod_vec_init(nres);

    for (i = 0; i < fmpz_mat_nrows(mat); i++)
    {
        for (j = 0; j < fmpz_mat_ncols(mat); j++)
        {
            for (k = 0; k < nres; k++)
                r[k] = nmod_mat_entry(residues[k], i, j);
            fmpz_multi_CRT_ui(fmpz_mat_entry(mat, i, j), r, comb, temp, sign);
        }
    }

    _nmod_vec_clear(r);
}

void
fmpz_mat_multi_CRT_ui(fmpz_mat_t mat, nmod_mat_t * const residues,
    slong nres, int sign)
{
    fmpz_comb_t comb;
    fmpz_comb_temp_t temp;
    mp_ptr primes;
    slong i;

    primes = _nmod_vec_init(nres);
    for (i = 0; i < nres; i++)
        primes[i] = residues[i]->mod.n;

    fmpz_comb_init(comb, primes, nres);
    fmpz_comb_temp_init(temp, comb);

    fmpz_mat_multi_CRT_ui_precomp(mat, residues, nres, comb, temp, sign);

    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(temp);
    _nmod_vec_clear(primes);
}
