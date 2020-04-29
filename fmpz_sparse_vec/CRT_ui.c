/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_sparse_vec.h"
#include "fmpz_sparse_vec.h"

void
fmpz_sparse_vec_CRT_ui(fmpz_sparse_vec_t w, const fmpz_sparse_vec_t u, const fmpz_t m1, 
                       const nmod_sparse_vec_t v, nmod_t m2, mp_limb_t m1i_m2, int sign)
{
    slong unnz = u->nnz, vnnz = v->nnz, wnnz, k;
    fmpz_sparse_entry_struct *ue, *we;
    nmod_sparse_entry_struct *ve;
    fmpz_t m1m2, zero;
    fmpz_init(m1m2);
    fmpz_init_set_ui(zero, UWORD(0));
    fmpz_mul_ui(m1m2, m1, m2.n);

    wnnz = _fmpz_sparse_vec_union_nnz_nmod (u, v);
    _fmpz_sparse_vec_resize(w, wnnz);
    ue = u->entries + unnz, ve = v->entries + vnnz, we = w->entries + wnnz;
    while ((k = _fmpz_sparse_vector_merge_descend_nmod (&we, &ue, &ve, u, v)) >= 0)
    {
        _fmpz_CRT_ui_precomp(we->val, (k==1)?zero:ue->val, m1, (k==0)?UWORD(0):ve->val, 
                             m2.n, m2.ninv, m1m2, m1i_m2, sign);
    }
    fmpz_clear(zero);
    fmpz_clear(m1m2);
}

