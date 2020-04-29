/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_vec.h"

/*         for (j = 0; j < fmpz_vec_ncols(vec); j++)
        {
            
            for (k = 0; k < nres; k++)
                nmod_vec_entry(residues[k], i, j) = r[k];
        }
 */
void
fmpz_sparse_vec_multi_mod_ui_precomp(nmod_sparse_vec_struct * residues, slong nres, const fmpz_sparse_vec_t v, 
                                     const fmpz_comb_t comb, fmpz_comb_temp_t temp)
{
    slong i, j;
    nmod_sparse_entry_struct *e;
    nmod_sparse_vec_struct *vm;
    mp_limb_t *r;

    r = flint_malloc(nres * sizeof(*r));

    for (j = 0; j < nres; ++j)
        vm = &residues[j], vm->entries = realloc(vm->entries, v->nnz * sizeof(*vm->entries)), vm->nnz = 0;

    for (i = 0; i < v->nnz; i++)
    {
        fmpz_multi_mod_ui(r, v->entries[i].val, comb, temp);
        for (j = 0; j < nres; ++j) 
            if (r[j] != UWORD(0))
                e = &residues[j].entries[residues[j].nnz++], e->ind = v->entries[i].ind, e->val = r[j];
    }

    for (j = 0; j < nres; ++j)
        vm = &residues[j], vm->entries = realloc(vm->entries, vm->nnz * sizeof(*vm->entries));

    flint_free(r);
}

void
fmpz_sparse_vec_multi_mod_ui(nmod_sparse_vec_struct * residues, mp_srcptr primes, slong nres, const fmpz_sparse_vec_t v)
{
    fmpz_comb_t comb;
    fmpz_comb_temp_t temp;

    fmpz_comb_init(comb, primes, nres);
    fmpz_comb_temp_init(temp, comb);

    fmpz_sparse_vec_multi_mod_ui_precomp(residues, nres, v, comb, temp);

    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(temp);
}
