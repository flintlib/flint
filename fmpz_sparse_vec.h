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

#ifndef FMPZ_SPARSE_VEC_H
#define FMPZ_SPARSE_VEC_H

#ifdef FMPZ_SPARSE_VEC_INLINES_C
#define FMPZ_SPARSE_VEC_INLINE FLINT_DLL
#else
#define FMPZ_SPARSE_VEC_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <string.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "nmod_sparse_vec.h"
#include "fmpz_vec.h"
#include "fmpz.h"
#include "thread_support.h"

#ifdef __cplusplus
  extern "C" {
#endif

typedef struct {
    slong ind;
    fmpz_t val;
} fmpz_sparse_entry_struct;

typedef fmpz_sparse_entry_struct fmpz_sparse_entry_t[1];

typedef struct {
    fmpz_sparse_entry_struct *entries;
    slong nnz;
} fmpz_sparse_vec_struct;

typedef fmpz_sparse_vec_struct fmpz_sparse_vec_t[1];

FMPZ_SPARSE_VEC_INLINE 
int fmpz_sparse_entry_cmp(const void *va, const void *vb)
{
    const fmpz_sparse_entry_struct *a = va;
    const fmpz_sparse_entry_struct *b = vb;
    if (a->ind < b->ind) return -1;
    if (b->ind < a->ind) return 1;
    return 0;
}

/* Memory management */
FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_init(fmpz_sparse_vec_t vec) 
{
    memset(vec, 0, sizeof(*vec));
}
FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_clear(fmpz_sparse_vec_t vec) 
{
    slong i;
    for (i = 0; i < vec->nnz; ++i)
        fmpz_clear(vec->entries[i].val);
    flint_free(vec->entries);
    memset(vec, 0, sizeof(*vec));
}

FMPZ_SPARSE_VEC_INLINE
void _fmpz_sparse_vec_resize(fmpz_sparse_vec_t vec, slong nnz) 
{
    slong i;
    if (nnz == 0) fmpz_sparse_vec_clear(vec);
    else if (nnz != vec->nnz)
    {
        for (i = nnz; i < vec->nnz; ++i)
            fmpz_clear (vec->entries[i].val);
        vec->entries = flint_realloc(vec->entries, nnz*sizeof(*vec->entries));
        for (i = vec->nnz; i < nnz; ++i)
            fmpz_init (vec->entries[i].val);
    }
    vec->nnz = nnz;
}

FMPZ_SPARSE_VEC_INLINE 
slong _fmpz_sparse_vec_support(slong **supp, fmpz_sparse_vec_t vec)
{
    slong i;
    *supp = flint_malloc(vec->nnz*sizeof(**supp));
    for (i = 0; i < vec->nnz; ++i) (*supp)[i] = vec->entries[i].ind;
    return vec->nnz;
}
FMPZ_SPARSE_VEC_INLINE 
void fmpz_sparse_vec_swap(fmpz_sparse_vec_t vec1, fmpz_sparse_vec_t vec2) 
{
    fmpz_sparse_vec_t tmp;
    *tmp = *vec1, *vec1 = *vec2, *vec2 = *tmp;
}

FLINT_DLL
slong fmpz_sparse_vec_max_bits(const fmpz_sparse_vec_t v);

/* Vector indexing */
FLINT_DLL 
fmpz_t * fmpz_sparse_vec_at(const fmpz_sparse_vec_t vec, slong i);

/* One-time instantiation */
FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_zero(fmpz_sparse_vec_t vec) 
{
    fmpz_sparse_vec_clear(vec);
}

FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_one(fmpz_sparse_vec_t vec, slong ind) 
{
    _fmpz_sparse_vec_resize(vec, 1);
    vec->entries[0].ind = ind;
    fmpz_one(vec->entries[0].val);
}

FLINT_DLL
void fmpz_sparse_vec_set(fmpz_sparse_vec_t vec, const fmpz_sparse_vec_t src, slong ioff);

FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_set_entry(fmpz_sparse_vec_t v, slong ind, const fmpz_t val)
{
    fmpz_t * oval = fmpz_sparse_vec_at(v, ind);
    if (oval == NULL)
    {
        _fmpz_sparse_vec_resize(v, v->nnz+1);
        v->entries[v->nnz-1].ind = ind;
        fmpz_set(v->entries[v->nnz-1].val, val);
        if (v->nnz >= 2 && v->entries[v->nnz-2].ind > ind)
            qsort(v->entries, v->nnz, sizeof(*v->entries), fmpz_sparse_entry_cmp);
    }
    else fmpz_set(*oval, val);
}

FLINT_DLL
void fmpz_sparse_vec_from_entries(fmpz_sparse_vec_t vec, slong * inds, fmpz * vals, slong nnz);

/* Vector comparison */
FMPZ_SPARSE_VEC_INLINE
int fmpz_sparse_vec_is_zero(const fmpz_sparse_vec_t vec) 
{
    return vec->nnz == 0;
}

FLINT_DLL
int fmpz_sparse_vec_equal(const fmpz_sparse_vec_t vec1, const fmpz_sparse_vec_t vec2, slong ioff);

/* Convert from/to dense vector */
FLINT_DLL
void fmpz_sparse_vec_from_dense(fmpz_sparse_vec_t vec, const fmpz *src, slong len);

FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_to_dense(fmpz *vec, const fmpz_sparse_vec_t src, slong len)
{
    slong i;
    _fmpz_vec_zero(vec, len);
    for (i = 0; i < src->nnz; ++i) 
        if (src->entries[i].ind < len) 
            fmpz_set(&vec[src->entries[i].ind], src->entries[i].val);
}

/* To/from modular vector(s) */
FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_get_nmod_sparse_vec(nmod_sparse_vec_t dst, const fmpz_sparse_vec_t src, const nmod_t mod)
{
    slong i, ind;
    if (src->nnz == 0)
        nmod_sparse_vec_zero(dst);
    else
    {
        dst->entries = flint_realloc(dst->entries, src->nnz*sizeof(*dst->entries));
        for (i = ind = 0; i < src->nnz; ++i) 
        {
            dst->entries[ind].ind = src->entries[i].ind;
            dst->entries[ind].val = fmpz_fdiv_ui(src->entries[i].val, mod.n);
            if (dst->entries[ind].val != UWORD(0)) ++ind;
        }
        if (ind == 0)
        {
            flint_free(dst->entries);
            dst->entries = NULL;
        }
        else dst->entries = flint_realloc(dst->entries, ind*sizeof(*dst->entries));
        dst->nnz = ind;
    }    
}

FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_set_nmod_sparse_vec_unsigned(fmpz_sparse_vec_t dst, const nmod_sparse_vec_t src)
{
    slong i;
    _fmpz_sparse_vec_resize(dst, src->nnz);
    for (i = 0; i < src->nnz; ++i) 
    {
        dst->entries[i].ind = src->entries[i].ind;
        fmpz_set_ui(dst->entries[i].val, src->entries[i].val);
    }
}

FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_set_nmod_sparse_vec(fmpz_sparse_vec_t dst, const nmod_sparse_vec_t src, const nmod_t mod)
{
    slong i;
    _fmpz_sparse_vec_resize(dst, src->nnz);
    for (i = 0; i < src->nnz; ++i) 
    {
        dst->entries[i].ind = src->entries[i].ind;
        fmpz_set_ui_smod(dst->entries[i].val, src->entries[i].val, mod.n);
    }
}

FLINT_DLL
void fmpz_sparse_vec_multi_mod_ui_precomp(nmod_sparse_vec_struct * residues, slong nres, const fmpz_sparse_vec_t v, 
                                          const fmpz_comb_t comb, fmpz_comb_temp_t temp);

FLINT_DLL
void fmpz_sparse_vec_multi_mod_ui(nmod_sparse_vec_struct * residues, mp_srcptr primes, slong nres, const fmpz_sparse_vec_t v);

FLINT_DLL
void fmpz_sparse_vec_CRT_ui(fmpz_sparse_vec_t w, const fmpz_sparse_vec_t u, const fmpz_t m1, 
                            const nmod_sparse_vec_t v, nmod_t m2, mp_limb_t m1i_m2, int sign);

FLINT_DLL
void fmpz_sparse_vec_multi_CRT_ui_precomp(fmpz_sparse_vec_t v, nmod_sparse_vec_struct const * residues, slong nres,
                                          const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign);

FLINT_DLL
void fmpz_sparse_vec_multi_CRT_ui(fmpz_sparse_vec_t v, nmod_sparse_vec_struct * const residues, mp_srcptr primes, slong nres, int sign);

/* Windows and concatenation */
FLINT_DLL
void fmpz_sparse_vec_window_init(fmpz_sparse_vec_t window, const fmpz_sparse_vec_t vec, slong i1, slong i2);


FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_window_clear(fmpz_sparse_vec_t window) 
{
    memset(window, 0, sizeof(*window));
}

FLINT_DLL
void fmpz_sparse_vec_concat(fmpz_sparse_vec_t res, const fmpz_sparse_vec_t vec1,  const fmpz_sparse_vec_t vec2, slong len1);

FLINT_DLL
void fmpz_sparse_vec_split(fmpz_sparse_vec_t res1, fmpz_sparse_vec_t res2, const fmpz_sparse_vec_t vec, slong ind);

/* Vector permutation */
FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_permute_inds(fmpz_sparse_vec_t vec, slong *P) 
{
    slong i;
    for (i = 0; i < vec->nnz; ++i) vec->entries[i].ind = P[vec->entries[i].ind];
}

/* Random vector generation */
FLINT_DLL 
void fmpz_sparse_vec_randtest(fmpz_sparse_vec_t vec, flint_rand_t state, slong nnz, slong len, flint_bitcnt_t bits);

FLINT_DLL 
void fmpz_sparse_vec_randtest_unsigned(fmpz_sparse_vec_t vec, flint_rand_t state, slong nnz, slong len, flint_bitcnt_t bits);

/* Vector display */
FLINT_DLL
void fmpz_sparse_vec_print_pretty(const fmpz_sparse_vec_t vec, slong ioff, slong maxi);

/* Vector operations */
FMPZ_SPARSE_VEC_INLINE
void fmpz_sparse_vec_neg(fmpz_sparse_vec_t v, const fmpz_sparse_vec_t u)
{
    slong i;
    fmpz_sparse_vec_set(v, u, 0);
    for (i = 0; i < v->nnz; ++i) fmpz_neg(v->entries[i].val, v->entries[i].val);
}

FLINT_DLL
void fmpz_sparse_vec_scalar_mul_fmpz(fmpz_sparse_vec_t v, const fmpz_sparse_vec_t u, const fmpz_t c);

FLINT_DLL
void fmpz_sparse_vec_scalar_divexact_fmpz(fmpz_sparse_vec_t v, const fmpz_sparse_vec_t u, const fmpz_t c);

FLINT_DLL
void fmpz_sparse_vec_scalar_mod_fmpz(fmpz_sparse_vec_t v, const fmpz_sparse_vec_t u, const fmpz_t mod);

FLINT_DLL
void fmpz_sparse_vec_scalar_mods_fmpz(fmpz_sparse_vec_t v, const fmpz_sparse_vec_t u, const fmpz_t mod);

/* Utility macros used by binary vector operations */
/* Compute total number of indices between two sparse vectors */
FMPZ_SPARSE_VEC_INLINE
slong _fmpz_sparse_vec_union_nnz(const fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v)
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

/* Utility macros used by binary vector operations */
/* Compute total number of indices between two sparse vectors */
FMPZ_SPARSE_VEC_INLINE
slong _fmpz_sparse_vec_union_nnz_nmod(const fmpz_sparse_vec_t u, const nmod_sparse_vec_t v)
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
FMPZ_SPARSE_VEC_INLINE
slong _fmpz_sparse_vector_merge_descend(fmpz_sparse_entry_struct **we, 
                                        fmpz_sparse_entry_struct **ue, 
                                        fmpz_sparse_entry_struct **ve, 
                                        const fmpz_sparse_vec_t u, 
                                        const fmpz_sparse_vec_t v)
{
    slong uind = (*ue==u->entries) ? -1 : (*ue-1)->ind;
    slong vind = (*ve==v->entries) ? -1 : (*ve-1)->ind;
    if (uind == -1 && vind == -1) return -1;
    if (uind == vind) {--*ue, --*ve, --*we; (*we)->ind = uind; return 2;}
    if (uind < vind) {--*ve, --*we; (*we)->ind = vind; return 1;}
    --*ue, --*we; (*we)->ind = uind; return 0;
}

/* Iterate through u and v in descending order, assigning sorted indices to w */
/* Returns -1 if u and v are both exhausted, 
    2 if we->ind = ue->ind == ve->ind
    1 if we->ind = ve->ind > ue->ind (or u exhausted), 
    0 if we->ind = ue->ind > ve->ind (or v exhausted). */
FMPZ_SPARSE_VEC_INLINE
slong _fmpz_sparse_vector_merge_descend_nmod(fmpz_sparse_entry_struct **we, 
                                        fmpz_sparse_entry_struct **ue, 
                                        nmod_sparse_entry_struct **ve, 
                                        const fmpz_sparse_vec_t u, 
                                        const nmod_sparse_vec_t v)
{
    slong uind = (*ue==u->entries) ? -1 : (*ue-1)->ind;
    slong vind = (*ve==v->entries) ? -1 : (*ve-1)->ind;
    if (uind == -1 && vind == -1) return -1;
    if (uind == vind) {--*ue, --*ve, --*we; (*we)->ind = uind; return 2;}
    if (uind < vind) {--*ve, --*we; (*we)->ind = vind; return 1;}
    --*ue, --*we; (*we)->ind = uind; return 0;
}

/* Like resize, but removes entries from the front of the vector */
FMPZ_SPARSE_VEC_INLINE
void _fmpz_sparse_vector_shift_left (fmpz_sparse_vec_t v, slong amt)
{
    slong i;
    if (amt == v->nnz) fmpz_sparse_vec_clear(v);
    else if (amt > 0)
    {
        for (i = amt; i < v->nnz; ++i)
        {
            v->entries[i-amt].ind = v->entries[i].ind;
            fmpz_set(v->entries[i - amt].val, v->entries[i].val);
        }
        _fmpz_sparse_vec_resize(v, v->nnz - amt);
    }
}

FLINT_DLL
void fmpz_sparse_vec_add(fmpz_sparse_vec_t w, const fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v);

FLINT_DLL
void fmpz_sparse_vec_sub(fmpz_sparse_vec_t w, const fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v);

FLINT_DLL 
void fmpz_sparse_vec_scalar_addmul_fmpz(fmpz_sparse_vec_t w, const fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v, const fmpz_t c);

FLINT_DLL
void fmpz_sparse_vec_scalar_submul_fmpz(fmpz_sparse_vec_t w, const fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v, const fmpz_t c);

FLINT_DLL
void fmpz_sparse_vec_dot(fmpz_t ret, const fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v);

FLINT_DLL
void fmpz_sparse_vec_dot_dense(fmpz_t ret, const fmpz_sparse_vec_t u, const fmpz *v);

/* Gaussian elimination operations */
FLINT_DLL
void fmpz_sparse_vec_gauss_elim_col(fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v, slong col);

FLINT_DLL
void fmpz_sparse_vec_gauss_elim(fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v);

FLINT_DLL
void fmpz_sparse_vec_gauss_elim_ext(fmpz_sparse_vec_t u, fmpz_sparse_vec_t v);

#ifdef __cplusplus
}
#endif

#endif

