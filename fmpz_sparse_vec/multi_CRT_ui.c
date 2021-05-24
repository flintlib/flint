/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_vec.h"

void
fmpz_sparse_vec_multi_CRT_ui_precomp(fmpz_sparse_vec_t v, nmod_sparse_vec_struct const * residues, slong nres,
                                     const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign)
{
    slong j, pos, *rpos, max_nnz, max_ind;
    mp_limb_t *r;
    fmpz_sparse_entry_struct *e;

    for (j = max_nnz = max_ind = 0; j < nres; ++j)
    {
        max_nnz = FLINT_MAX(max_nnz, residues[j].nnz);
        if (residues[j].nnz != 0) max_ind = FLINT_MAX(max_ind, residues[j].entries[residues[j].nnz - 1].ind);
    }
    _fmpz_sparse_vec_resize(v, max_nnz); /* May change later */
    if (max_nnz == 0) return;

    rpos = flint_calloc(nres, sizeof(*rpos));
    r = flint_malloc(nres * sizeof (*r));
    for (pos = 0; ; ++pos)
    {
        if (pos == max_nnz) max_nnz *= 2, _fmpz_sparse_vec_resize(v, max_nnz);
        e = &v->entries[pos];

        /* Get next minimal index */
        e->ind = max_ind + 1;
        for (j = 0; j < nres; ++j)
        {
            e->ind = (rpos[j] != residues[j].nnz && residues[j].entries[rpos[j]].ind < e->ind) ? residues[j].entries[rpos[j]].ind : e->ind;
        }
        if (e->ind == max_ind + 1) break;
        for (j = 0; j < nres; ++j)
        {
            r[j] = (rpos[j] != residues[j].nnz && residues[j].entries[rpos[j]].ind == e->ind) ? residues[j].entries[rpos[j]++].val : 0;
        }
        fmpz_multi_CRT_ui(e->val, r, comb, temp, sign);
    }
    _fmpz_sparse_vec_resize(v, pos);
    flint_free(rpos);
    flint_free(r);
}

void
fmpz_sparse_vec_multi_CRT_ui(fmpz_sparse_vec_t v, nmod_sparse_vec_struct * const residues, mp_srcptr primes, slong nres, int sign)
{
    fmpz_comb_t comb;
    fmpz_comb_temp_t temp;

    fmpz_comb_init(comb, primes, nres);
    fmpz_comb_temp_init(temp, comb);

    fmpz_sparse_vec_multi_CRT_ui_precomp(v, residues, nres, comb, temp, sign);

    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(temp);
}
