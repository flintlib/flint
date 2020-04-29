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
#include "fmpz_sparse_vec.h"

void fmpz_sparse_vec_scalar_addmul_fmpz(fmpz_sparse_vec_t w, const fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v, const fmpz_t c)
{
    slong unnz = u->nnz, vnnz = v->nnz, wnnz, k;
    fmpz_t tmp;
    fmpz_sparse_entry_struct *ue, *ve, *we;

    /* Check for simpler operations first */
    if (fmpz_is_zero(c) || vnnz == 0) fmpz_sparse_vec_set(w, u, 0);
    else if (fmpz_is_one(c)) fmpz_sparse_vec_add(w, u, v);
    else if (fmpz_equal_si(c, WORD(-1))) fmpz_sparse_vec_sub(w, u, v);
    else if (unnz == 0) fmpz_sparse_vec_scalar_mul_fmpz(w, v, c);
    else
    {
        fmpz_init(tmp);
        wnnz = _fmpz_sparse_vec_union_nnz (u, v);
        _fmpz_sparse_vec_resize(w, wnnz);
        ue = u->entries + unnz, ve = v->entries + vnnz, we = w->entries + wnnz;
        while ((k = _fmpz_sparse_vector_merge_descend (&we, &ue, &ve, u, v)) >= 0)
        {
            switch(k)
            {
                case 0: fmpz_set(we->val, ue->val); break;
                case 1: fmpz_mul(we->val, ve->val, c); break;
                default: fmpz_mul(tmp, ve->val, c);
                         fmpz_add(we->val, ue->val, tmp);
                         if (fmpz_is_zero(we->val)) we++;
            }
        }
        _fmpz_sparse_vector_shift_left (w, we - w->entries);
        fmpz_clear(tmp);
    }
}
