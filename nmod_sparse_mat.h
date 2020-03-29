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

#ifndef NMOD_SPARSE_MAT_H
#define NMOD_SPARSE_MAT_H

#ifdef NMOD_SPARSE_MAT_INLINES_C
#define NMOD_SPARSE_MAT_INLINE FLINT_DLL
#else
#define NMOD_SPARSE_MAT_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_sparse_vec.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "thread_support.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* A sparse matrix is a list of sparse vectors */
typedef struct
{
    nmod_sparse_vec_struct *rows;
    slong r;
    slong c;
    slong c_off;
    nmod_t mod;
}
nmod_sparse_mat_struct;

typedef nmod_sparse_mat_struct nmod_sparse_mat_t[1];

NMOD_SPARSE_MAT_INLINE
void _nmod_sparse_mat_set_mod(nmod_sparse_mat_t mat, mp_limb_t n)
{
    mat->mod.n = n;
    count_leading_zeros(mat->mod.norm, n);
    invert_limb(mat->mod.ninv, n << mat->mod.norm);
}

/* Memory management */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_init(nmod_sparse_mat_t mat, slong rows, slong cols, nmod_t mod) 
{
    mat->rows = flint_calloc(rows, sizeof(*mat->rows));
    mat->r = rows;
    mat->c = cols;
    mat->c_off = 0;
    mat->mod = mod;
}
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_clear(nmod_sparse_mat_t mat) 
{
    slong i;
    for(i=0; i<mat->r; ++i) nmod_sparse_vec_clear(&mat->rows[i]);
    flint_free(mat->rows);
    memset(mat, 0, sizeof(*mat));
}
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_swap(nmod_sparse_mat_t mat1, nmod_sparse_mat_t mat2) 
{
    nmod_sparse_mat_t tmp;
    *tmp = *mat1; *mat1 = *mat2; *mat2 = *tmp;
}

/* One-time instantiation */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_zero(nmod_sparse_mat_t mat) 
{
    slong i;
    for(i=0; i<mat->r; ++i) nmod_sparse_vec_zero(&mat->rows[i]);
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_one(nmod_sparse_mat_t mat) 
{
    slong i;
    for(i=0; i<mat->r; ++i) nmod_sparse_vec_one(&mat->rows[i], i);
}

FLINT_DLL
void nmod_sparse_mat_set(nmod_sparse_mat_t mat, const nmod_sparse_mat_t src);

FLINT_DLL
void nmod_sparse_mat_from_entries(nmod_sparse_mat_t mat, slong * rows, slong * cols, mp_limb_t * vals, slong nnz);

/* Convert from/to dense matrix */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_from_dense(nmod_sparse_mat_t mat, const nmod_mat_t src)
{
    slong i;
    for(i=0; i<src->r; ++i) nmod_sparse_vec_from_dense(&mat->rows[i], src->rows[i], src->c);
}
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_to_dense(nmod_mat_t mat, const nmod_sparse_mat_t src) {
    slong i;
    for(i=0; i<src->r; ++i) nmod_sparse_vec_to_dense(mat->rows[i], &src->rows[i], mat->c);
}

/* Windows and concatenation */
FLINT_DLL
void nmod_sparse_mat_window_init(nmod_sparse_mat_t window, const nmod_sparse_mat_t mat, slong r1, slong c1, slong r2, slong c2);

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_window_clear(nmod_sparse_mat_t window) 
{
    flint_free(window->rows);
    memset(window, 0, sizeof(*window));
}
<<<<<<< HEAD

/* res->r must equal mat1->r and mat2->r */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_concat_horizontal(nmod_sparse_mat_t res,
                                    const nmod_sparse_mat_t mat1,  const nmod_sparse_mat_t mat2) 
{
    slong i;
    res->c = mat1->c + mat2->c;
    for(i=0; i<res->r; ++i)
        nmod_sparse_vec_concat(&res->rows[i], &mat1->rows[i], &mat2->rows[i], mat1->c);
}
/* res->r must equal mat1->r + mat2->r */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_concat_vertical(nmod_sparse_mat_t res, const nmod_sparse_mat_t mat1,  const nmod_sparse_mat_t mat2) 
{
    slong i;
    res->c = FLINT_MAX(mat1->c, mat2->c);
    for(i=0; i<mat1->r; ++i)
        nmod_sparse_vec_set(&res->rows[i], &mat1->rows[i]);
    for(i=mat1->r; i<res->r; ++i)
        nmod_sparse_vec_set(&res->rows[i], &mat2->rows[i-mat1->r]);
}

=======

/* res->r must equal mat1->r and mat2->r */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_concat_horizontal(nmod_sparse_mat_t res,
                                    const nmod_sparse_mat_t mat1,  const nmod_sparse_mat_t mat2) 
{
    slong i;
    res->c = mat1->c + mat2->c;
    for(i=0; i<res->r; ++i)
        nmod_sparse_vec_concat(&res->rows[i], &mat1->rows[i], &mat2->rows[i], mat1->c);
}
/* res->r must equal mat1->r + mat2->r */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_concat_vertical(nmod_sparse_mat_t res, const nmod_sparse_mat_t mat1,  const nmod_sparse_mat_t mat2) 
{
    slong i;
    res->c = FLINT_MAX(mat1->c, mat2->c);
    for(i=0; i<mat1->r; ++i)
        nmod_sparse_vec_set(&res->rows[i], &mat1->rows[i]);
    for(i=mat1->r; i<res->r; ++i)
        nmod_sparse_vec_set(&res->rows[i], &mat2->rows[i-mat1->r]);
}

>>>>>>> Added sparse vector class to nmod, changed sparse matrix class to use it for underlying, added (untested) LU decomposition
/* Matrix permutation */
NMOD_SPARSE_VEC_INLINE
void nmod_sparse_mat_permute_cols(nmod_sparse_mat_t mat, slong *Q) 
{
    slong i;
    for(i=0; i<mat->r; ++i) nmod_sparse_vec_permute_inds(&mat->rows[i], Q);
}

/* Random matrix generation */
FLINT_DLL void nmod_sparse_mat_randtest(nmod_sparse_mat_t mat, flint_rand_t state, slong min_nnz, slong max_nnz);
/*
FLINT_DLL void nmod_sparse_mat_randfull(nmod_sparse_mat_t mat, flint_rand_t state);
FLINT_DLL int nmod_sparse_mat_randpermdiag(nmod_sparse_mat_t mat, flint_rand_t state,
                 mp_srcptr diag, slong n);
FLINT_DLL void nmod_sparse_mat_randrank(nmod_sparse_mat_t, flint_rand_t state, slong rank);
FLINT_DLL void nmod_sparse_mat_randops(nmod_sparse_mat_t mat, slong count, flint_rand_t state);
FLINT_DLL void nmod_sparse_mat_randtril(nmod_sparse_mat_t mat, flint_rand_t state, int unit);
FLINT_DLL void nmod_sparse_mat_randtriu(nmod_sparse_mat_t mat, flint_rand_t state, int unit);
 */

FLINT_DLL void nmod_sparse_mat_print_pretty(const nmod_sparse_mat_t mat);

NMOD_SPARSE_MAT_INLINE
int nmod_sparse_mat_equal(const nmod_sparse_mat_t mat1, const nmod_sparse_mat_t mat2) 
{
    slong i;
    if(mat1->r != mat2->r) return 0;
    for(i=0; i<mat1->r; ++i)
        if(nmod_sparse_vec_equal(&mat1->rows[i], &mat2->rows[i], mat1->c_off-mat2->c_off)==0) return 0;
    return 1;
}

NMOD_SPARSE_MAT_INLINE
int nmod_sparse_mat_is_zero(const nmod_sparse_mat_t mat) 
{
    slong i;
    for(i=0; i<mat->r; ++i) 
        if(!nmod_sparse_vec_is_zero(&mat->rows[i])) return 0;
    return 1;
}

NMOD_SPARSE_MAT_INLINE
int nmod_sparse_mat_is_square(const nmod_sparse_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Must have A->r == B->c and A->c == B->r */
FLINT_DLL void nmod_sparse_mat_transpose(nmod_sparse_mat_t B, const nmod_sparse_mat_t A);

/* Arithmetic */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_neg(nmod_sparse_mat_t B, const nmod_sparse_mat_t A) 
{
    slong i;
    nmod_sparse_mat_set(B, A);
    for(i=0; i<B->r; ++i) nmod_sparse_vec_neg(&B->rows[i], &B->rows[i], B->mod);
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_scalar_mul(nmod_sparse_mat_t B, const nmod_sparse_mat_t A, mp_limb_t c) 
{
    if(c==UWORD(0)) nmod_sparse_mat_zero(B);
    else {
        slong i;
        nmod_sparse_mat_set(B, A);
        for(i=0; i<B->r; ++i) nmod_sparse_vec_scalar_mul(&B->rows[i], &B->rows[i], c, B->mod);    
    }
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_scalar_mul_fmpz(nmod_sparse_mat_t res, const nmod_sparse_mat_t M, const fmpz_t c)
{
    fmpz_t d;
    fmpz_init(d);
    fmpz_mod_ui(d, c, res->mod.n);
    nmod_sparse_mat_scalar_mul(res, M, fmpz_get_ui(d));
    fmpz_clear(d);
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_addmul(nmod_sparse_mat_t C, const nmod_sparse_mat_t A, const nmod_sparse_mat_t B, mp_limb_t c) 
{
    slong i;
    for(i=0; i<C->r; ++i) nmod_sparse_vec_scalar_addmul(&C->rows[i], &A->rows[i], &B->rows[i], c, C->mod);
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_add(nmod_sparse_mat_t C, const nmod_sparse_mat_t A, const nmod_sparse_mat_t B) 
{
    nmod_sparse_mat_addmul(C, A, B, UWORD(1));
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_sub(nmod_sparse_mat_t C, const nmod_sparse_mat_t A, const nmod_sparse_mat_t B) 
{
    nmod_sparse_mat_addmul(C, A, B, C->mod.n-UWORD(1));
}

/* Matrix-vector and matrix-matrix multipliciation */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_mul_vec(mp_ptr y, const nmod_sparse_mat_t A, const mp_ptr x) 
{
    slong i;
    for(i=0; i<A->r; ++i) y[i] = nmod_sparse_vec_dot_dense(&A->rows[i], x, A->mod);
}
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_mul_mat(nmod_mat_t Y, const nmod_sparse_mat_t A, const nmod_mat_t X) 
{
    slong i, j;
    for(i=0; i<A->r; ++i) {
        for(j=0; j<A->rows[i].nnz; ++j) {
            nmod_sparse_entry_struct *e = &A->rows[i].entries[j];
            _nmod_vec_scalar_addmul_nmod(Y->rows[i], X->rows[e->ind], X->c, e->val, Y->mod);
        }
    }
}

/* Permutations */
/* FLINT_DLL void nmod_sparse_mat_swap_rows(nmod_sparse_mat_t mat, slong * perm, slong r, slong s);
FLINT_DLL void nmod_sparse_mat_invert_rows(nmod_sparse_mat_t mat, slong * perm);
FLINT_DLL void nmod_sparse_mat_swap_cols(nmod_sparse_mat_t mat, slong * perm, slong r, slong s);
FLINT_DLL void nmod_sparse_mat_invert_cols(nmod_sparse_mat_t mat, slong * perm);
FLINT_DLL void nmod_sparse_mat_apply_permutation(nmod_sparse_mat_t A, slong * P, slong n);
 */

/* Decomposition */
void nmod_sparse_mat_lu(slong *P, slong *Q, nmod_sparse_mat_t L, nmod_sparse_mat_t U, const nmod_sparse_mat_t A);

/* Nonsingular solving */
int nmod_sparse_mat_solve_lanczos(mp_ptr x, const nmod_sparse_mat_t A, const mp_ptr b, flint_rand_t state);
int nmod_sparse_mat_solve_lu(mp_ptr x, const nmod_sparse_mat_t A, const mp_ptr b);

/* Nullspace */
<<<<<<< HEAD
/* NMOD_SPARSE_MAT_INLINE
slong nmod_sparse_mat_nullspace(nmod_mat_t X, const nmod_sparse_mat_t A);
 */
=======
NMOD_SPARSE_MAT_INLINE
slong nmod_sparse_mat_nullspace(nmod_mat_t X, const nmod_sparse_mat_t A);

>>>>>>> Added sparse vector class to nmod, changed sparse matrix class to use it for underlying, added (untested) LU decomposition
#ifdef __cplusplus
}
#endif

#endif

