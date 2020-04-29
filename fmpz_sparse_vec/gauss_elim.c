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

/* Reduce u by v */
void fmpz_sparse_vec_gauss_elim_col(fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v, slong col)
{
    fmpz_t q;
    fmpz_t *uc = fmpz_sparse_vec_at(u, col);
    fmpz_t *vc = fmpz_sparse_vec_at(v, col);
    if (uc == NULL || vc == NULL) return;
    
    fmpz_init(q);
    fmpz_fdiv_q(q, *uc, *vc);
    fmpz_sparse_vec_scalar_submul_fmpz(u, u, v, q);
    fmpz_clear(q);
}

/* Reduce u by v */
void fmpz_sparse_vec_gauss_elim(fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v)
{
    fmpz_t q, *uc;
    fmpz_sparse_entry_struct *lu = u->entries, *lv = v->entries;
    if (u->nnz == 0 || v->nnz == 0 || lu->ind > lv->ind) return;
    fmpz_init(q);
    if (lu->ind == lv->ind) fmpz_fdiv_q(q, lu->val, lv->val);
    else if ((uc = fmpz_sparse_vec_at(u, lv->ind))) fmpz_fdiv_q(q, *uc, lv->val);
    fmpz_sparse_vec_scalar_submul_fmpz(u, u, v, q);
    fmpz_clear(q);
}

/* Apply unimodular transformation to (u,v) to minimize both vectors lexicographically */
void fmpz_sparse_vec_gauss_elim_ext(fmpz_sparse_vec_t u, fmpz_sparse_vec_t v)
{
    slong vnnz = v->nnz, unnz = u->nnz, nnz, k, pc, i;
    fmpz_sparse_entry_struct *lu = u->entries;
    fmpz_sparse_entry_struct *lv = v->entries;
    fmpz_t g, vv, vu, uv, uu, a, b;
    fmpz_sparse_entry_struct *ue, *ve, *nue, *nve;
    if (u->nnz == 0 || v->nnz == 0) return;
    if (lu->ind != lv->ind) {fmpz_sparse_vec_gauss_elim(u, v); return;}
    pc = lu->ind;
    if (fmpz_cmpabs(lu->val, lv->val) < 0) /* Pre-apply transform [[0, -1], [1, 0]] */
    {
        fmpz_sparse_vec_swap(u, v);
        fmpz_sparse_vec_neg(u, u);
        lu = u->entries, lv = v->entries;
        vnnz = v->nnz, unnz = u->nnz;
    }
    if (fmpz_sgn(lv->val) < 0) /* Pre-apply transform [[-1, 0], [0, -1]] */
    {
        fmpz_sparse_vec_neg(v, v);
        fmpz_sparse_vec_neg(u, u);
    }

    /* Check for trivial cases */
    if (fmpz_divisible(lu->val, lv->val)) {fmpz_sparse_vec_gauss_elim(u, v); return;}

    fmpz_init(g);
    fmpz_init(vv);
    fmpz_init(vu);
    fmpz_init(uv);
    fmpz_init(uu);
    fmpz_init(a);
    fmpz_init(b);

    /* Construct transformation */
    fmpz_xgcd(g, vv, vu, lv->val, lu->val);
    if (fmpz_sgn(g) < 0)
    {
        fmpz_neg(vv, vv);
        fmpz_neg(vu, vu);
    }
    fmpz_divexact(uu, lv->val, g);
    fmpz_divexact(uv, lu->val, g); 
    fmpz_neg(uv, uv); /* [[uu uv] [vu vv]] is unimodular */

    /* Reallocate vectors */
    nnz = _fmpz_sparse_vec_union_nnz (u, v);
    _fmpz_sparse_vec_resize(u, nnz);
    _fmpz_sparse_vec_resize(v, nnz);
    ue = u->entries + unnz, ve = v->entries + vnnz; /* Old locations */
    nue = u->entries + nnz, nve = v->entries + nnz; /* New locations */
    while ((k = _fmpz_sparse_vector_merge_descend (&nue, &ue, &ve, u, v)) >= 0)
    {
        nve--;
        switch(k)
        {
            case 0: nve->ind = ue->ind; fmpz_mul(nve->val, ue->val, vu); fmpz_mul(nue->val, ue->val, uu); break;
            case 1: nve->ind = ve->ind; fmpz_mul(nue->val, ve->val, uv); fmpz_mul(nve->val, ve->val, vv); break;
            default: nve->ind = ve->ind; 
            fmpz_set(a, ve->val);
            fmpz_set(b, ue->val); 
            fmpz_mul(nve->val, a, vv); fmpz_addmul(nve->val, b, vu);
            fmpz_mul(nue->val, b, uu); fmpz_addmul(nue->val, a, uv);
        }
        if (fmpz_is_zero(nue->val)) nue++;
        if (fmpz_is_zero(nve->val)) nve++;
    }
    _fmpz_sparse_vector_shift_left (u, nue - u->entries);
    _fmpz_sparse_vector_shift_left (v, nve - v->entries);
    fmpz_clear(g);
    fmpz_clear(vv);
    fmpz_clear(vu);
    fmpz_clear(uv);
    fmpz_clear(uu);
    fmpz_clear(a);
    fmpz_clear(b);
}
