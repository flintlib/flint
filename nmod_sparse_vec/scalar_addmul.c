/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_vec.h"

void nmod_sparse_vec_scalar_addmul(nmod_sparse_vec_t w, const nmod_sparse_vec_t u, const nmod_sparse_vec_t v, mp_limb_t c, nmod_t mod)
{
    if (c == UWORD(0) || (u->nnz == 0 && v->nnz == 0)) {nmod_sparse_vec_zero(w); return;}
    if (v->nnz == 0) {nmod_sparse_vec_set(w, u, 0); return;}
    if (u->nnz == 0) {nmod_sparse_vec_scalar_mul(w, v, c, mod); return;}

    nmod_sparse_vec_t u2; u2->nnz = u->nnz;
    if (u==w)
    {
        /* Handle inplace operations */
        u2->entries = flint_malloc(u->nnz*sizeof(*u2->entries));
        memcpy(u2->entries, u->entries, u->nnz*sizeof(*u2->entries));
    } 
    else u2->entries = u->entries;

    w->entries = flint_realloc(w->entries, (u->nnz+v->nnz)*sizeof(*w->entries));
    memset(w->entries, 0, (u->nnz+v->nnz)*sizeof(*w->entries));
    slong i = 0, j = 0, k = 0;
    nmod_sparse_entry_struct *ue = u2->entries;
    nmod_sparse_entry_struct *ve = v->entries;
    nmod_sparse_entry_struct *we = w->entries;

    /* Interleave u and v */
    while (i < u->nnz && j < v->nnz)
    {
        we[k].ind = FLINT_MIN(ue[i].ind, ve[j].ind);
        if (ue[i].ind == we[k].ind) we[k].val = ue[i].val, ++i;
        if (ve[j].ind == we[k].ind) we[k].val = nmod_addmul(we[k].val, ve[j].val, c, mod), ++j;
        if (we[k].val != UWORD(0)) ++k;
    }
    while (i < u->nnz) we[k].ind = ue[i].ind, we[k].val = ue[i].val, ++i, ++k;
    while (j < v->nnz) we[k].ind = ve[j].ind, we[k].val = nmod_mul(ve[j].val, c, mod), ++j, ++k;
    
    if (u==w) flint_free(u2->entries);
    if (k == 0) nmod_sparse_vec_clear(w);
    else w->entries = realloc(w->entries, k*sizeof(*w->entries));
    w->nnz = k;
}
