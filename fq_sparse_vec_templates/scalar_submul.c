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

void TEMPLATE(T, TEMPLATE(sparse_vec_scalar_submul, T))(TEMPLATE(T, sparse_vec_t) w, const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, sparse_vec_t) v, const TEMPLATE(T, t) c, 
                const TEMPLATE(T, ctx_t) ctx)
{
    slong i, nnz, unnz, vnnz;
    TEMPLATE(T, t) tmp;
    TEMPLATE(T, sparse_entry_struct) *ue, *ve, *we;
    if (u->nnz == 0) 
    {
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T)) (w, v, c, ctx);
        TEMPLATE(T, sparse_vec_neg) (w, w, ctx);
    }
    else if (v->nnz == 0 || TEMPLATE(T, is_zero) (c, ctx)) TEMPLATE(T, sparse_vec_set) (w, u, 0, ctx);
    else 
    {
        TEMPLATE(T, init) (tmp, ctx);
        unnz = u->nnz; vnnz = v->nnz; /* Store in case either equals w (gets resized) */
        nnz = _TEMPLATE(T, sparse_vec_union_nnz) (u, v, ctx);
        _TEMPLATE(T, sparse_vec_resize) (w, nnz, ctx);
        ue = u->entries + unnz - 1, ve = v->entries + vnnz - 1, we = w->entries + w->nnz - 1;
        while (ue >= u->entries && ve >= v->entries)
        {
            we->ind = FLINT_MAX(ue->ind, ve->ind);
            if (ue->ind == ve->ind) 
            {
                TEMPLATE(T, mul) (tmp, (ve--)->val, c, ctx);
                TEMPLATE(T, sub) (we->val, tmp, (ue--)->val, ctx);
                if (!TEMPLATE(T, is_zero) (we->val, ctx)) we--;
            }
            else if (ue->ind == we->ind) TEMPLATE(T, set) ((we--)->val, (ue--)->val, ctx);
            else if (ve->ind == we->ind) 
            {
                TEMPLATE(T, mul) (we->val, (ve--)->val, c, ctx);
                TEMPLATE(T, neg) (we->val, we->val, c, ctx); we--;
            }
        }
        while (ue >= u->entries) we->ind = ue->ind, TEMPLATE(T, set) ((we--)->val, (ue--)->val, ctx);
        while (ve >= v->entries) 
        {
            we->ind = ve->ind;
            TEMPLATE(T, mul) (we->val, (ve--)->val, c, ctx);
            TEMPLATE(T, neg) (we->val, we->val, c, ctx);
            we--;
        }

        nnz = w->nnz - (we + 1 - w->entries);
        if (nnz == 0) TEMPLATE(T, sparse_vec_clear) (w, ctx);
        else if (nnz < w->nnz)
        {
            for (i = 0; i < w->nnz - nnz; ++i) TEMPLATE(T, clear) (w->entries[i].val, ctx);
            memmove(w->entries, we + 1, nnz*sizeof(*w->entries));
            w->entries = flint_realloc(w->entries, nnz*sizeof(*w->entries));
            w->nnz = nnz;
        }
        TEMPLATE(T, clear) (tmp, ctx);
    }
}
#endif

