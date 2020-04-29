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
#include "flint.h"
#include "fmpz_sparse_vec.h"

void fmpz_sparse_vec_add(fmpz_sparse_vec_t w, const fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v)
{
    slong unnz = u->nnz, vnnz = v->nnz, wnnz, k;
    fmpz_sparse_entry_struct *ue, *ve, *we;

    /* Check for simpler operations first */
    if (vnnz == 0) fmpz_sparse_vec_set(w, u, 0);
    else if (unnz == 0) fmpz_sparse_vec_set(w, v, 0);
    else
    {
        wnnz = _fmpz_sparse_vec_union_nnz (u, v);
        _fmpz_sparse_vec_resize(w, wnnz);
        ue = u->entries + unnz, ve = v->entries + vnnz, we = w->entries + wnnz;
        while ((k = _fmpz_sparse_vector_merge_descend (&we, &ue, &ve, u, v)) >= 0)
        {
            switch(k)
            {
                case 0: fmpz_set(we->val, ue->val); break;
                case 1: fmpz_set(we->val, ve->val); break;
                default: fmpz_add(we->val, ue->val, ve->val);
                         if (fmpz_is_zero(we->val)) we++;
            }
        }
        _fmpz_sparse_vector_shift_left (w, we - w->entries);
    }
}
