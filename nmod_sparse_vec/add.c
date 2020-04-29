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

void nmod_sparse_vec_add(nmod_sparse_vec_t w, const nmod_sparse_vec_t u, const nmod_sparse_vec_t v, nmod_t mod)
{
    slong unnz = u->nnz, vnnz = v->nnz, wnnz, k;
    nmod_sparse_entry_struct *ue, *ve, *we;

    /* Check for simpler operations first */
    if (vnnz == 0) nmod_sparse_vec_set(w, u, 0);
    else if (unnz == 0) nmod_sparse_vec_set(w, v, 0);
    else
    {
        wnnz = _nmod_sparse_vec_union_nnz (u, v);
        w->entries = flint_realloc(w->entries, wnnz*sizeof(*w->entries));
        w->nnz = wnnz;
        ue = u->entries + unnz, ve = v->entries + vnnz, we = w->entries + wnnz;
        while ((k = _nmod_sparse_vector_merge_descend (&we, &ue, &ve, u, v)) >= 0)
        {
            switch(k)
            {
                case 0: we->val = ue->val; break;
                case 1: we->val = ve->val; break;
                default: we->val = nmod_add(ue->val, ve->val, mod);
                         if (we->val == UWORD(0)) we++;
            }
        }
        _nmod_sparse_vector_shift_left (w, we - w->entries);
    }
}
