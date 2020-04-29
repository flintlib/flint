/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010,2011 Fredrik Johansson
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include <string.h>
#include "flint.h"
#include "longlong.h"
#include "templates.h"
#include "ulong_extras.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct {
    slong ind;
    TEMPLATE(T, t) val;
} TEMPLATE(T, sparse_entry_struct);

typedef TEMPLATE(T, sparse_entry_struct) TEMPLATE(T, sparse_entry_t)[1];

typedef struct {
    TEMPLATE(T, sparse_entry_struct) *entries;
    slong nnz;
} TEMPLATE(T, sparse_vec_struct);

typedef TEMPLATE(T, sparse_vec_struct) TEMPLATE(T, sparse_vec_t)[1];

FQ_SPARSE_VEC_TEMPLATES_INLINE
int TEMPLATE(T, sparse_entry_cmp)(const void *va, const void *vb)
{
    const TEMPLATE(T, sparse_entry_struct) *a = va;
    const TEMPLATE(T, sparse_entry_struct) *b = vb;
    if (a->ind < b->ind) return -1;
    if (b->ind < a->ind) return 1;
    return 0;
}

/* Memory management */
FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_init)(TEMPLATE(T, sparse_vec_t) vec, 
                                    const TEMPLATE(T, ctx_t) ctx) 
{
    memset(vec, 0, sizeof(*vec));
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_clear)(TEMPLATE(T, sparse_vec_t) vec, 
                                    const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < vec->nnz; ++i)
        TEMPLATE(T, clear) (vec->entries[i].val, ctx);
    flint_free(vec->entries);
    memset(vec, 0, sizeof(*vec));
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void _TEMPLATE(T, sparse_vec_resize)(TEMPLATE(T, sparse_vec_t) vec, slong nnz,
                                    const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    if (nnz == 0) TEMPLATE(T, sparse_vec_clear) (vec, ctx);
    else if (nnz != vec->nnz)
    {
        for (i = nnz; i < vec->nnz; ++i)
            TEMPLATE(T, clear) (vec->entries[i].val, ctx);
        vec->entries = flint_realloc(vec->entries, nnz*sizeof(*vec->entries));
        for (i = vec->nnz; i < nnz; ++i)
            TEMPLATE(T, init) (vec->entries[i].val, ctx);
    }
    vec->nnz = nnz;
}

FQ_SPARSE_VEC_TEMPLATES_INLINE 
void TEMPLATE(T, sparse_vec_swap)(TEMPLATE(T, sparse_vec_t) vec1, TEMPLATE(T, sparse_vec_t) vec2, 
                                    const TEMPLATE(T, ctx_t) ctx) 
{
    TEMPLATE(T, sparse_vec_t) tmp;
    *tmp = *vec1, *vec1 = *vec2, *vec2 = *tmp;
}

/* Vector indexing */
FQ_SPARSE_VEC_TEMPLATES_INLINE
TEMPLATE(T, t) *TEMPLATE(T, sparse_vec_at)(TEMPLATE(T, sparse_vec_t) vec, slong i,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong j;
    for (j = 0; j < vec->nnz; ++j)
        if (vec->entries[j].ind==i) 
            return &vec->entries[j].val;
    return NULL;
}

/* One-time instantiation */
FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_zero)(TEMPLATE(T, sparse_vec_t) vec, 
                                    const TEMPLATE(T, ctx_t) ctx) 
{
    TEMPLATE(T, sparse_vec_clear)(vec, ctx);
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_one)(TEMPLATE(T, sparse_vec_t) vec, slong ind, 
                                    const TEMPLATE(T, ctx_t) ctx) 
{
    _TEMPLATE(T, sparse_vec_resize) (vec, 1, ctx);
    vec->entries[0].ind = ind;
    TEMPLATE(T, one) (vec->entries[0].val, ctx);
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_set)(TEMPLATE(T, sparse_vec_t) dst, const TEMPLATE(T, sparse_vec_t) src, slong ioff,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    if (dst == src) return;
    _TEMPLATE(T, sparse_vec_resize) (dst, src->nnz, ctx);
    for (i = 0; i < dst->nnz; ++i)
    {
        dst->entries[i].ind = src->entries[i].ind - ioff;
        TEMPLATE(T, set) (dst->entries[i].val, src->entries[i].val, ctx);
    }
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_set_entry)(TEMPLATE(T, sparse_vec_t) v, slong ind, const TEMPLATE(T, t) val, 
                                    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) *oval;
    if (TEMPLATE(T, is_zero) (val, ctx)) return;
    oval = TEMPLATE(T, sparse_vec_at) (v, ind, ctx);
    if (oval == NULL)
    {
        _TEMPLATE(T, sparse_vec_resize) (v, v->nnz + 1, ctx);
        TEMPLATE(T, set) (v->entries[v->nnz-1].val, val, ctx);
        v->entries[v->nnz-1].ind = ind;
        if (v->nnz >= 2 && ind < v->entries[v->nnz-2].ind)
            qsort(v->entries, v->nnz, sizeof(*v->entries), TEMPLATE(T, sparse_entry_cmp));
    }
    else TEMPLATE(T, set) (*oval, val, ctx);
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_from_entries)(TEMPLATE(T, sparse_vec_t) vec, slong * inds, TEMPLATE(T, struct) * vals, slong nnz,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    _TEMPLATE(T, sparse_vec_resize) (vec, nnz, ctx);
    for (i = 0; i < nnz; ++i)
    {
        vec->entries[i].ind = inds[i];
        TEMPLATE(T, set) (vec->entries[i].val, &vals[i], ctx);
    }
}

/* Vector comparison */
FQ_SPARSE_VEC_TEMPLATES_INLINE
int TEMPLATE(T, sparse_vec_is_zero)(const TEMPLATE(T, sparse_vec_t) vec, 
                                    const TEMPLATE(T, ctx_t) ctx) 
{
    return vec->nnz == 0;
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
int TEMPLATE(T, sparse_vec_equal)(const TEMPLATE(T, sparse_vec_t) vec1, const TEMPLATE(T, sparse_vec_t) vec2, slong ioff,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    if (vec1->nnz != vec2->nnz) return 0;
    for (i = 0; i < vec1->nnz; ++i)
    {
        if ((vec1->entries[i].ind != vec2->entries[i].ind + ioff) || 
           !TEMPLATE(T, equal)(vec1->entries[i].val, vec2->entries[i].val, ctx)) return 0;
    }
    return 1;
}

/* Convert from/to dense vector */
FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_from_dense)(TEMPLATE(T, sparse_vec_t) dst, TEMPLATE(T, struct) *src, slong len,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, nnz = 0;
    TEMPLATE(T, sparse_entry_struct) *e;
    for (i = 0; i < len; ++i)
        if (!TEMPLATE(T, is_zero) (&src[i], ctx)) ++nnz;
    _TEMPLATE(T, sparse_vec_resize) (dst, nnz, ctx);
    e = dst->entries;
    for (i = 0; i < len; ++i)
        if (!TEMPLATE(T, is_zero) (&src[i], ctx))
            e->ind = i, TEMPLATE(T, set) (e->val, &src[i], ctx), ++e;
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_to_dense)(TEMPLATE(T, struct) *dst, const TEMPLATE(T, sparse_vec_t) src, slong len, 
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    _TEMPLATE(T, vec_zero)(dst, len, ctx);
    for (i = 0; i < src->nnz; ++i) TEMPLATE(T, set)(&dst[src->entries[i].ind], src->entries[i].val, ctx);
}

/* Windows and concatenation */
FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_window_init)(TEMPLATE(T, sparse_vec_t) window, const TEMPLATE(T, sparse_vec_t) vec, slong i1, slong i2,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong start, end;
    for (start = 0; start < vec->nnz && vec->entries[start].ind < i1; ++start);
    for (end = vec->nnz; end > 0 && vec->entries[end-1].ind >= i2; --end);
    window->entries = vec->entries + start;
    window->nnz = end - start;
}


FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_window_clear)(TEMPLATE(T, sparse_vec_t) window, 
                                    const TEMPLATE(T, ctx_t) ctx) 
{
    memset(window, 0, sizeof(*window));
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_concat)(TEMPLATE(T, sparse_vec_t) res, const TEMPLATE(T, sparse_vec_t) vec1,  const TEMPLATE(T, sparse_vec_t) vec2, slong len1,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, nnz = vec1->nnz+vec2->nnz;
    _TEMPLATE(T, sparse_vec_resize) (res, nnz, ctx);
    for (i = 0; i < nnz; ++i)
    {
        TEMPLATE(T, sparse_entry_struct) *e = (i < vec1->nnz) ? &vec1->entries[i] : &vec2->entries[i-vec1->nnz];
        res->entries[i].ind = e->ind + ((i < vec1->nnz) ? 0 : len1);
        TEMPLATE(T, set) (res->entries[i].val, e->val, ctx);
    }
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_split)(TEMPLATE(T, sparse_vec_t) res1, TEMPLATE(T, sparse_vec_t) res2, const TEMPLATE(T, sparse_vec_t) vec, slong ind,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, nnz1;
    TEMPLATE(T, sparse_entry_struct) *e1, *e2, *e;
    for (nnz1 = 0; nnz1 < vec->nnz; ++nnz1) if (vec->entries[nnz1].ind >= ind) break;

    _TEMPLATE(T, sparse_vec_resize) (res1, nnz1, ctx);
    _TEMPLATE(T, sparse_vec_resize) (res2, vec->nnz - nnz1, ctx);
    e1 = res1->entries, e2 = res2->entries;
    for (i = 0; i < vec->nnz; ++i)
    {
        e = (i < nnz1) ? e1++ : e2++;
        e->ind = vec->entries[i].ind - ((i < nnz1) ?  0 : ind);
        TEMPLATE(T, set) (e->val, vec->entries[i].val, ctx);
    }
}

/* Vector permutation */
FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_permute_inds)(TEMPLATE(T, sparse_vec_t) vec, slong *P, 
                                    const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < vec->nnz; ++i) vec->entries[i].ind = P[vec->entries[i].ind];
    qsort(vec->entries, vec->nnz, sizeof(*vec->entries), TEMPLATE(T, sparse_entry_cmp));
}

/* Random vector generation */
FLINT_DLL 
void TEMPLATE(T, sparse_vec_randtest)(TEMPLATE(T, sparse_vec_t) vec, flint_rand_t state, slong nnz, slong len,
                                    const TEMPLATE(T, ctx_t) ctx);

/* Vector display */
FLINT_DLL
void TEMPLATE(T, sparse_vec_print_pretty)(const TEMPLATE(T, sparse_vec_t) vec, slong ioff, slong maxi,
                                    const TEMPLATE(T, ctx_t) ctx);

/* Vector operations */
FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_neg)(TEMPLATE(T, sparse_vec_t) v, const TEMPLATE(T, sparse_vec_t) u,  
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    TEMPLATE(T, sparse_vec_set)(v, u, 0, ctx);
    for (i = 0; i < v->nnz; ++i) TEMPLATE(T, neg) (v->entries[i].val, v->entries[i].val, ctx);
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T))(TEMPLATE(T, sparse_vec_t) v, const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, t) c,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    if (TEMPLATE(T, is_zero) (c, ctx)) {TEMPLATE(T, sparse_vec_zero)(v, ctx); return;}
    TEMPLATE(T, sparse_vec_set)(v, u, 0, ctx);
    if (!TEMPLATE(T, is_one)(c, ctx))
        for (i = 0; i < v->nnz; ++i) TEMPLATE(T, mul)(v->entries[i].val, v->entries[i].val, c, ctx);
}

/* Utility macros used by binary vector operations */
/* Compute total number of indices between two sparse vectors */
FQ_SPARSE_VEC_TEMPLATES_INLINE
slong _TEMPLATE(T, sparse_vec_union_nnz)(const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, sparse_vec_t) v,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, nnz = 0;
    for (i = j = 0; i < u->nnz && j < v->nnz; ++nnz)
    {
        if (u->entries[i].ind == v->entries[j].ind) ++i, ++j;
        else if (u->entries[i].ind < v->entries[j].ind) ++i;
        else if (u->entries[i].ind > v->entries[j].ind) ++j;
    }
    nnz += u->nnz - i + v->nnz - j;
    return nnz;
}

/* Iterate through u and v in descending order, assigning sorted indices to w */
/* Returns -1 if u and v are both exhausted, 
    2 if we->ind = ue->ind == ve->ind
    1 if we->ind = ve->ind > ue->ind (or u exhausted), 
    0 if we->ind = ue->ind > ve->ind (or v exhausted). */
FQ_SPARSE_VEC_TEMPLATES_INLINE
slong _TEMPLATE(T, sparse_vector_merge_descend) (TEMPLATE(T, sparse_entry_struct) **we, 
                                                 TEMPLATE(T, sparse_entry_struct) **ue, 
                                                 TEMPLATE(T, sparse_entry_struct) **ve, 
                                                 const TEMPLATE(T, sparse_vec_t) u, 
                                                 const TEMPLATE(T, sparse_vec_t) v)
{
    slong uind = (*ue==u->entries) ? -1 : (*ue-1)->ind;
    slong vind = (*ve==v->entries) ? -1 : (*ve-1)->ind;
    if (uind == -1 && vind == -1) return -1;
    if (uind == vind) {--*ue, --*ve, --*we; (*we)->ind = uind; return 2;}
    if (uind < vind) {--*ve, --*we; (*we)->ind = vind; return 1;}
    --*ue, --*we; (*we)->ind = uind; return 0;
}

/* Like resize, but removes entries from the front of the vector */
FQ_SPARSE_VEC_TEMPLATES_INLINE
void _TEMPLATE(T, sparse_vector_shift_left) (TEMPLATE(T, sparse_vec_t) v, slong amt, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    if (amt == v->nnz) TEMPLATE(T, sparse_vec_clear) (v, ctx);
    else if (amt > 0)
    {
        v->nnz -= amt;
        for (i = 0; i < amt; ++i) TEMPLATE(T, clear) (v->entries[i].val, ctx);
        memmove(v->entries, v->entries + amt, v->nnz*sizeof(*v->entries));
        v->entries = flint_realloc(v->entries, v->nnz*sizeof(*v->entries));
    }
}

FLINT_DLL 
void TEMPLATE(T, sparse_vec_add)(TEMPLATE(T, sparse_vec_t) w, const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, sparse_vec_t) v,  
                                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL 
void TEMPLATE(T, sparse_vec_sub)(TEMPLATE(T, sparse_vec_t) w, const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, sparse_vec_t) v,  
                                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL 
void TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T))(TEMPLATE(T, sparse_vec_t) w, const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, sparse_vec_t) v, const TEMPLATE(T, t) c,
                                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL 
void TEMPLATE(T, TEMPLATE(sparse_vec_scalar_submul, T))(TEMPLATE(T, sparse_vec_t) w, const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, sparse_vec_t) v, const TEMPLATE(T, t) c,
                                    const TEMPLATE(T, ctx_t) ctx);

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_dot)(TEMPLATE(T, t) ret, const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, sparse_vec_t) v,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;
    TEMPLATE(T, t) tmp;
    TEMPLATE(T, init) (tmp, ctx);
    TEMPLATE(T, zero) (ret, ctx);
    for (i = j = 0; i < u->nnz && j < v->nnz; )
    {
        if (u->entries[i].ind == v->entries[j].ind)
        {
            TEMPLATE(T, mul) (tmp, u->entries[i].val, v->entries[j].val, ctx);
            TEMPLATE(T, add) (ret, ret, tmp, ctx), ++i, ++j;
        }
        else if (u->entries[i].ind < v->entries[j].ind) ++i;
        else if (u->entries[i].ind > v->entries[j].ind) ++j;
    }
    TEMPLATE(T, clear) (tmp, ctx);
}

FQ_SPARSE_VEC_TEMPLATES_INLINE
void TEMPLATE(T, sparse_vec_dot_dense)(TEMPLATE(T, t) ret, const TEMPLATE(T, sparse_vec_t) u, const TEMPLATE(T, struct) *v,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    TEMPLATE(T, t) tmp;
    TEMPLATE(T, init) (tmp, ctx);
    TEMPLATE(T, zero) (ret, ctx);
    for (i = 0; i < u->nnz; ++i) 
    {
        TEMPLATE(T, mul) (tmp, u->entries[i].val, &v[u->entries[i].ind],  ctx);
        TEMPLATE(T, add) (ret, ret, tmp, ctx);
    }
    TEMPLATE(T, clear) (tmp, ctx);
}

#ifdef __cplusplus
}
#endif

#endif

