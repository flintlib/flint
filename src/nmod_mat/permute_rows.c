/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
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
    mp_limb_t ** mat_tmp = (mp_limb_t **) flint_malloc(mat->r * sizeof(mp_limb_t *));

    /* perm_store[i] <- perm_store[perm_act[i]] */
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->r);

    /* rows[i] <- rows[perm_act[i]]  */
    for (i = 0; i < mat->r; i++)
        mat_tmp[i] = mat->rows[perm_act[i]];
    for (i = 0; i < mat->r; i++)
        mat->rows[i] = mat_tmp[i];

    flint_free(mat_tmp);
}
