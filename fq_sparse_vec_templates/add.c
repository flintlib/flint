/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
#ifdef T

#include <string.h>
#include "templates.h"

void TEMPLATE(T, sparse_vec_add)(TEMPLATE(T, sparse_vec_t) w, const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, sparse_vec_t) v,
                const TEMPLATE(T, ctx_t) ctx)
{
    slong unnz = u->nnz, vnnz = v->nnz, wnnz, k;
    TEMPLATE(T, sparse_entry_struct) *ue, *ve, *we;

    if (vnnz == 0) {TEMPLATE(T, sparse_vec_set) (w, u, 0, ctx); return;}
    if (unnz == 0) {TEMPLATE(T, sparse_vec_set) (w, v, 0, ctx); return;}
    wnnz = _TEMPLATE(T, sparse_vec_union_nnz) (u, v, ctx);
    _TEMPLATE(T, sparse_vec_resize) (w, wnnz, ctx);
    ue = u->entries + unnz, ve = v->entries + vnnz, we = w->entries + wnnz;
    while ((k = _TEMPLATE(T, sparse_vector_merge_descend) (&we, &ue, &ve, u, v)) >= 0)
    {
        switch(k)
        {
        case 0: TEMPLATE(T, set)(we->val, ue->val, ctx); break;
        case 1: TEMPLATE(T, set)(we->val, ve->val, ctx); break;
        default: TEMPLATE(T, add)(we->val, ue->val, ve->val, ctx); 
                 if (TEMPLATE(T, is_zero) (we->val, ctx)) we++;
        }
    }
    _TEMPLATE(T, sparse_vector_shift_left) (w, we - w->entries, ctx);
}
#endif

