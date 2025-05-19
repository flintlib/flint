/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

/** Permute rows of a matrix `mat` according to `perm_act`, and propagate the
 * action on `perm_store`.
 * That is, performs for each appropriate index `i`, the operations
 * `perm_store[i] <- perm_store[perm_act[i]]`
 * `rows[i] <- rows[perm_act[i]]` */
void
nmod_mat_permute_rows(nmod_mat_t mat, const slong * perm_act, slong * perm_store)
{
    slong i;
    ulong * mat_tmp = (ulong *) flint_malloc(mat->r * mat->c * sizeof(ulong));

    /* perm_store[i] <- perm_store[perm_act[i]] */
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->r);

    /* rows[i] <- rows[perm_act[i]]  */
    for (i = 0; i < mat->r; i++)
        _nmod_vec_set(mat_tmp + i * mat->c, nmod_mat_entry_ptr(mat, perm_act[i], 0), mat->c);

    for (i = 0; i < mat->r; i++)
        _nmod_vec_set(nmod_mat_entry_ptr(mat, i, 0), mat_tmp + i * mat->c, mat->c);

    flint_free(mat_tmp);
}
