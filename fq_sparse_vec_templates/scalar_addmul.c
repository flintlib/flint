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

void TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T))(TEMPLATE(T, sparse_vec_t) w, const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, sparse_vec_t) v, const TEMPLATE(T, t) c, 
                const TEMPLATE(T, ctx_t) ctx)
{
    slong unnz = u->nnz, vnnz = v->nnz, wnnz, k;
    TEMPLATE(T, sparse_entry_struct) *ue, *ve, *we;
    TEMPLATE(T, t) tmp;

    /* Check for simpler operations first */
    if (vnnz == 0 || TEMPLATE(T, is_zero) (c, ctx)) {TEMPLATE(T, sparse_vec_set) (w, u, 0, ctx); return;}
    if (TEMPLATE(T, is_one) (c, ctx)) {TEMPLATE(T, sparse_vec_add) (w, u, v, ctx); return;}
    if (unnz == 0) {TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T)) (w, v, c, ctx); return;}
    TEMPLATE(T, init) (tmp, ctx);
    TEMPLATE(T, neg) (tmp, c, ctx);
    if (TEMPLATE(T, is_one) (tmp, ctx)) TEMPLATE(T, sparse_vec_sub) (w, u, v, ctx);
    else /* Now just do standard addmul */ 
    {
        wnnz = _TEMPLATE(T, sparse_vec_union_nnz) (u, v, ctx);
        _TEMPLATE(T, sparse_vec_resize) (w, wnnz, ctx);
        ue = u->entries + unnz, ve = v->entries + vnnz, we = w->entries + wnnz;
        while ((k = _TEMPLATE(T, sparse_vector_merge_descend) (&we, &ue, &ve, u, v)) >= 0)
        {
            switch(k)
            {
                case 0: TEMPLATE(T, set)(we->val, ue->val, ctx); break;
                case 1: TEMPLATE(T, mul)(we->val, ve->val, c, ctx); break;
                default: TEMPLATE(T, mul) (tmp, ve->val, c, ctx);
                         TEMPLATE(T, add) (we->val, ue->val, tmp, ctx);
                         if (TEMPLATE(T, is_zero) (we->val, ctx)) we++;
            }
        }
        _TEMPLATE(T, sparse_vector_shift_left) (w, we - w->entries, ctx);
    }
    
    TEMPLATE(T, clear) (tmp, ctx);
}
#endif

