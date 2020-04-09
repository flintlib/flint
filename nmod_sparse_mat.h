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

/* Memory management */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_init(nmod_sparse_mat_t M, slong rows, slong cols, nmod_t mod) 
{
    M->rows = flint_calloc(rows, sizeof(*M->rows));
    M->r = rows;
    M->c = cols;
    M->c_off = 0;
    M->mod = mod;
}
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_clear(nmod_sparse_mat_t M) 
{
    slong i;
    for (i = 0; i < M->r; ++i) nmod_sparse_vec_clear(&M->rows[i]);
    flint_free(M->rows);
    memset(M, 0, sizeof(*M));
}
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_swap(nmod_sparse_mat_t M1, nmod_sparse_mat_t M2) 
{
    nmod_sparse_mat_t tmp;
    *tmp = *M1; *M1 = *M2; *M2 = *tmp;
}

/* One-time instantiation */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_zero(nmod_sparse_mat_t M) 
{
    slong i;
    for (i = 0; i < M->r; ++i) nmod_sparse_vec_zero(&M->rows[i]);
}
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_one(nmod_sparse_mat_t M) 
{
    slong i;
    for (i = 0; i < M->r; ++i) nmod_sparse_vec_one(&M->rows[i], i);
}

/* One-time instantiation */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_set(nmod_sparse_mat_t M, const nmod_sparse_mat_t src) 
{
    slong i, rmax = FLINT_MIN(M->r, src->r);
    if (M==src || M->r == 0) return;
    for (i = 0; i < rmax; ++i) nmod_sparse_vec_set(&M->rows[i], &src->rows[i], src->c_off);
}

FLINT_DLL
void nmod_sparse_mat_from_entries(nmod_sparse_mat_t M, slong * rows, slong * cols, mp_limb_t * vals, slong nnz);

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_append_col(nmod_sparse_mat_t M, mp_srcptr v) 
{
    slong i;
    for (i = 0; i < M->r; ++i) nmod_sparse_vec_set_entry(&M->rows[i], M->c, v[i]);
    M->c += 1;
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_append_row(nmod_sparse_mat_t M, const nmod_sparse_vec_t v) 
{
    M->rows = realloc(M->rows, (M->r+1)*sizeof(*M->rows));
    memset(M->rows + M->r, 0, sizeof(*M->rows));
    nmod_sparse_vec_set(&M->rows[M->r], v, 0);
    M->r += 1;
}

/* Convert from/to dense matrix */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_from_dense(nmod_sparse_mat_t M, const nmod_mat_t dM)
{
    slong i, rmax = FLINT_MIN(M->r, dM->r);
    for (i = 0; i < rmax; ++i) nmod_sparse_vec_from_dense(&M->rows[i], dM->rows[i], dM->c);
}
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_to_dense(nmod_mat_t dM, const nmod_sparse_mat_t M)
{
    slong i, rmax = FLINT_MIN(M->r, dM->r);
    for (i = 0; i < rmax; ++i) nmod_sparse_vec_to_dense(dM->rows[i], &M->rows[i], dM->c);
}

/* Windows, concatenation, and splitting */
FLINT_DLL
void nmod_sparse_mat_window_init(nmod_sparse_mat_t W, const nmod_sparse_mat_t M, slong r1, slong c1, slong r2, slong c2);

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_window_clear(nmod_sparse_mat_t W) 
{
    flint_free(W->rows);
    memset(W, 0, sizeof(*W));
}


/* Combine M1 and M2 into block matrix B = [M1 M2] */
/* B->r must equal M1->r and M2->r */
<<<<<<< HEAD
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_concat_horizontal(nmod_sparse_mat_t B,
                                    const nmod_sparse_mat_t M1,  const nmod_sparse_mat_t M2) 
{
    slong i;
    B->c = M1->c + M2->c;
    for (i = 0; i < B->r; ++i)
        nmod_sparse_vec_concat(&B->rows[i], &M1->rows[i], &M2->rows[i], M1->c);
}
/* Combine M1 and M2 into block matrix B = [M1^t M1^t]^t */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_concat_vertical(nmod_sparse_mat_t B, const nmod_sparse_mat_t M1,  const nmod_sparse_mat_t M2) 
{
    slong i;
    B->c = FLINT_MAX(M1->c, M2->c);
    for (i = 0; i < M1->r; ++i)
        nmod_sparse_vec_set(&B->rows[i], &M1->rows[i], M1->c_off);
    for (i = M1->r; i < B->r; ++i)
        nmod_sparse_vec_set(&B->rows[i], &M2->rows[i-M1->r], M2->c_off);
}
<<<<<<< HEAD
<<<<<<< HEAD

/* Split block matrix B = [M1 M2] into submatrices M1 and M2 */
/* M1->r and M2->r must equal B->r */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_split_horizontal(nmod_sparse_mat_t M1, nmod_sparse_mat_t M2, const nmod_sparse_mat_t B, slong c)
{
    slong i;
    for (i = 0; i < B->r; ++i) nmod_sparse_vec_split(&M1->rows[i], &M2->rows[i], &B->rows[i], c);
}

/* Split block matix B = [M1^t M1^t]^t into submatrices M1 and M2 */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_split_vertical(nmod_sparse_mat_t M1, nmod_sparse_mat_t M2, const nmod_sparse_mat_t B, slong r)
{
    slong i;
    r = FLINT_MIN(r, B->r);
    for (i = 0; i < r; ++i) nmod_sparse_vec_set(&M1->rows[i], &B->rows[i], B->c_off);
    for (i = r; i < B->r; ++i) nmod_sparse_vec_set(&M2->rows[i-r], &B->rows[i], B->c_off);
}

<<<<<<< HEAD

=======
=======
=======
>>>>>>> Spacing and cuddling fixed


/* res->r must equal mat1->r and mat2->r */
=======
>>>>>>> Created sparse vector and matrix utilities for all the fq variants
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_concat_horizontal(nmod_sparse_mat_t B,
                                    const nmod_sparse_mat_t M1,  const nmod_sparse_mat_t M2) 
{
    slong i;
    B->c = M1->c + M2->c;
    for (i = 0; i < B->r; ++i)
        nmod_sparse_vec_concat(&B->rows[i], &M1->rows[i], &M2->rows[i], M1->c);
}
/* Combine M1 and M2 into block matrix B = [M1^t M1^t]^t */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_concat_vertical(nmod_sparse_mat_t B, const nmod_sparse_mat_t M1,  const nmod_sparse_mat_t M2) 
{
    slong i;
    B->c = FLINT_MAX(M1->c, M2->c);
    for (i = 0; i < M1->r; ++i)
        nmod_sparse_vec_set(&B->rows[i], &M1->rows[i], M1->c_off);
    for (i = M1->r; i < B->r; ++i)
        nmod_sparse_vec_set(&B->rows[i], &M2->rows[i-M1->r], M2->c_off);
}

/* Split block matrix B = [M1 M2] into submatrices M1 and M2 */
/* M1->r and M2->r must equal B->r */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_split_horizontal(nmod_sparse_mat_t M1, nmod_sparse_mat_t M2, const nmod_sparse_mat_t B, slong c)
{
    slong i;
    for (i = 0; i < B->r; ++i) nmod_sparse_vec_split(&M1->rows[i], &M2->rows[i], &B->rows[i], c);
}

/* Split block matix B = [M1^t M1^t]^t into submatrices M1 and M2 */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_split_vertical(nmod_sparse_mat_t M1, nmod_sparse_mat_t M2, const nmod_sparse_mat_t B, slong r)
{
    slong i;
    r = FLINT_MIN(r, B->r);
    for (i = 0; i < r; ++i) nmod_sparse_vec_set(&M1->rows[i], &B->rows[i], B->c_off);
    for (i = r; i < B->r; ++i) nmod_sparse_vec_set(&M2->rows[i-r], &B->rows[i], B->c_off);
}

<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> Added sparse vector class to nmod, changed sparse matrix class to use it for underlying, added (untested) LU decomposition
>>>>>>> Added sparse vector class to nmod, changed sparse matrix class to use it for underlying, added (untested) LU decomposition
=======
>>>>>>> Spacing and cuddling fixed
=======

>>>>>>> Now with additional utilities, more correct basic functions, and nullspace and inversion functions
/* Matrix permutation */
NMOD_SPARSE_VEC_INLINE
void nmod_sparse_mat_permute_cols(nmod_sparse_mat_t M, slong *Q) 
{
    slong i;
<<<<<<< HEAD
<<<<<<< HEAD
    for (i = 0; i < M->r; ++i) 
    {
        if (!M->rows[i].nnz) continue;
=======
    for (i = 0; i < M->r; ++i) {
        if(!M->rows[i].nnz) continue;
>>>>>>> Created sparse vector and matrix utilities for all the fq variants
        nmod_sparse_vec_permute_inds(&M->rows[i], Q);
        qsort(M->rows[i].entries, M->rows[i].nnz, sizeof(*M->rows[i].entries), nmod_sparse_entry_cmp);
    }
}

NMOD_SPARSE_VEC_INLINE
void nmod_sparse_mat_permute_rows(nmod_sparse_mat_t M, slong *P) 
{
    slong i;
    nmod_sparse_vec_struct *prows;
    prows = flint_calloc(M->r, sizeof(*prows));
    for (i = 0; i < M->r; ++i) prows[P[i]] = M->rows[i];
    memcpy(M->rows, prows, M->r*sizeof(*M->rows));
    flint_free(prows);
=======
    for (i = 0; i < mat->r; ++i) nmod_sparse_vec_permute_inds(&mat->rows[i], Q);
>>>>>>> Spacing and cuddling fixed
}

/* Random matrix generation */
FLINT_DLL void nmod_sparse_mat_randtest(nmod_sparse_mat_t M, flint_rand_t state, slong min_nnz, slong max_nnz);
/*
FLINT_DLL void nmod_sparse_mat_randfull(nmod_sparse_mat_t M, flint_rand_t state);
FLINT_DLL int nmod_sparse_mat_randpermdiag(nmod_sparse_mat_t M, flint_rand_t state,
                 mp_srcptr diag, slong n);
FLINT_DLL void nmod_sparse_mat_randrank(nmod_sparse_mat_t, flint_rand_t state, slong rank);
FLINT_DLL void nmod_sparse_mat_randops(nmod_sparse_mat_t M, slong count, flint_rand_t state);
FLINT_DLL void nmod_sparse_mat_randtril(nmod_sparse_mat_t M, flint_rand_t state, int unit);
FLINT_DLL void nmod_sparse_mat_randtriu(nmod_sparse_mat_t M, flint_rand_t state, int unit);
 */

FLINT_DLL void nmod_sparse_mat_print_pretty(const nmod_sparse_mat_t M);

NMOD_SPARSE_MAT_INLINE
int nmod_sparse_mat_equal(const nmod_sparse_mat_t M1, const nmod_sparse_mat_t M2) 
{
    slong i;
<<<<<<< HEAD
<<<<<<< HEAD
    if (M1->r != M2->r) return 0;
    for (i = 0; i < M1->r; ++i)
        if (nmod_sparse_vec_equal(&M1->rows[i], &M2->rows[i], M1->c_off-M2->c_off) == 0) return 0;
=======
    if (mat1->r != mat2->r) return 0;
    for (i = 0; i < mat1->r; ++i)
        if (nmod_sparse_vec_equal(&mat1->rows[i], &mat2->rows[i], mat1->c_off-mat2->c_off) == 0) return 0;
>>>>>>> Spacing and cuddling fixed
=======
    if (M1->r != M2->r) return 0;
    for (i = 0; i < M1->r; ++i)
        if (nmod_sparse_vec_equal(&M1->rows[i], &M2->rows[i], M1->c_off-M2->c_off) == 0) return 0;
>>>>>>> Created sparse vector and matrix utilities for all the fq variants
    return 1;
}

NMOD_SPARSE_MAT_INLINE
int nmod_sparse_mat_is_zero(const nmod_sparse_mat_t M) 
{
    slong i;
<<<<<<< HEAD
<<<<<<< HEAD
    for (i = 0; i < M->r; ++i) 
        if (!nmod_sparse_vec_is_zero(&M->rows[i])) return 0;
=======
    for (i = 0; i < mat->r; ++i) 
        if (!nmod_sparse_vec_is_zero(&mat->rows[i])) return 0;
>>>>>>> Spacing and cuddling fixed
    return 1;
}

=======
    for (i = 0; i < M->r; ++i) 
        if (!nmod_sparse_vec_is_zero(&M->rows[i])) return 0;
    return 1;
}

NMOD_SPARSE_MAT_INLINE
int nmod_sparse_mat_is_square(const nmod_sparse_mat_t M)
{
    return (M->r == M->c);
}

>>>>>>> Created sparse vector and matrix utilities for all the fq variants
/* Must have M->r == N->c and M->c == N->r */
FLINT_DLL void nmod_sparse_mat_transpose(nmod_sparse_mat_t N, const nmod_sparse_mat_t M);

/* Arithmetic */
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_neg(nmod_sparse_mat_t N, const nmod_sparse_mat_t M) 
{
    slong i;
<<<<<<< HEAD
<<<<<<< HEAD
    nmod_sparse_mat_set(N, M);
    for (i = 0; i < N->r; ++i) nmod_sparse_vec_neg(&N->rows[i], &N->rows[i], N->mod);
=======
    nmod_sparse_mat_set(B, A);
    for (i = 0; i < B->r; ++i) nmod_sparse_vec_neg(&B->rows[i], &B->rows[i], B->mod);
>>>>>>> Spacing and cuddling fixed
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_scalar_mul_nmod(nmod_sparse_mat_t N, const nmod_sparse_mat_t M, mp_limb_t c) 
{
<<<<<<< HEAD
    if (c == UWORD(0)) nmod_sparse_mat_zero(N);
    else {
        slong i;
        nmod_sparse_mat_set(N, M);
        for (i = 0; i < N->r; ++i) nmod_sparse_vec_scalar_mul_nmod(&N->rows[i], &N->rows[i], c, N->mod);    
=======
    if (c == UWORD(0)) nmod_sparse_mat_zero(B);
    else {
        slong i;
        nmod_sparse_mat_set(B, A);
        for (i = 0; i < B->r; ++i) nmod_sparse_vec_scalar_mul(&B->rows[i], &B->rows[i], c, B->mod);    
>>>>>>> Spacing and cuddling fixed
=======
    nmod_sparse_mat_set(N, M);
    for (i = 0; i < N->r; ++i) nmod_sparse_vec_neg(&N->rows[i], &N->rows[i], N->mod);
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_scalar_mul(nmod_sparse_mat_t N, const nmod_sparse_mat_t M, mp_limb_t c) 
{
    if (c == UWORD(0)) nmod_sparse_mat_zero(N);
    else {
        slong i;
        nmod_sparse_mat_set(N, M);
        for (i = 0; i < N->r; ++i) nmod_sparse_vec_scalar_mul(&N->rows[i], &N->rows[i], c, N->mod);    
>>>>>>> Created sparse vector and matrix utilities for all the fq variants
    }
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_scalar_mul_fmpz(nmod_sparse_mat_t N, const nmod_sparse_mat_t M, const fmpz_t c)
{
    fmpz_t d;
    fmpz_init(d);
    fmpz_mod_ui(d, c, N->mod.n);
<<<<<<< HEAD
    nmod_sparse_mat_scalar_mul_nmod(N, M, fmpz_get_ui(d));
=======
    nmod_sparse_mat_scalar_mul(N, M, fmpz_get_ui(d));
>>>>>>> Created sparse vector and matrix utilities for all the fq variants
    fmpz_clear(d);
}

NMOD_SPARSE_MAT_INLINE
<<<<<<< HEAD
void nmod_sparse_mat_add(nmod_sparse_mat_t O, const nmod_sparse_mat_t M, const nmod_sparse_mat_t N) 
{
    slong i;
    for (i = 0; i < O->r; ++i) nmod_sparse_vec_add(&O->rows[i], &M->rows[i], &N->rows[i], O->mod);
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_sub(nmod_sparse_mat_t O, const nmod_sparse_mat_t M, const nmod_sparse_mat_t N) 
{
    slong i;
    for (i = 0; i < O->r; ++i) nmod_sparse_vec_sub(&O->rows[i], &M->rows[i], &N->rows[i], O->mod);
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_scalar_addmul_nmod(nmod_sparse_mat_t O, const nmod_sparse_mat_t M, const nmod_sparse_mat_t N, mp_limb_t c) 
{
    slong i;
    for (i = 0; i < O->r; ++i) nmod_sparse_vec_scalar_addmul_nmod(&O->rows[i], &M->rows[i], &N->rows[i], c, O->mod);
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_scalar_submul_nmod(nmod_sparse_mat_t O, const nmod_sparse_mat_t M, const nmod_sparse_mat_t N, mp_limb_t c) 
{
    slong i;
    for (i = 0; i < O->r; ++i) nmod_sparse_vec_scalar_submul_nmod(&O->rows[i], &M->rows[i], &N->rows[i], c, O->mod);
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_addmul(nmod_sparse_mat_t C, const nmod_sparse_mat_t A, const nmod_sparse_mat_t B, mp_limb_t c) 
=======
void nmod_sparse_mat_addmul(nmod_sparse_mat_t O, const nmod_sparse_mat_t M, const nmod_sparse_mat_t N, mp_limb_t c) 
>>>>>>> Created sparse vector and matrix utilities for all the fq variants
{
    slong i;
    for (i = 0; i < O->r; ++i) nmod_sparse_vec_scalar_addmul(&O->rows[i], &M->rows[i], &N->rows[i], c, O->mod);
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_add(nmod_sparse_mat_t O, const nmod_sparse_mat_t M, const nmod_sparse_mat_t N) 
{
    nmod_sparse_mat_addmul(O, M, N, UWORD(1));
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_sub(nmod_sparse_mat_t O, const nmod_sparse_mat_t M, const nmod_sparse_mat_t N) 
{
    nmod_sparse_mat_addmul(O, M, N, O->mod.n-UWORD(1));
}

/* Matrix-vector and matrix-matrix multipliciation */
NMOD_SPARSE_MAT_INLINE
<<<<<<< HEAD
<<<<<<< HEAD
void nmod_sparse_mat_mul_vec(mp_ptr y, const nmod_sparse_mat_t M, mp_srcptr x) 
{
    slong i;
<<<<<<< HEAD
    for (i = 0; i < M->r; ++i) y[i] = nmod_sparse_vec_dot_dense(&M->rows[i], x, M->mod);
=======
    for (i = 0; i < A->r; ++i) y[i] = nmod_sparse_vec_dot_dense(&A->rows[i], x, A->mod);
>>>>>>> Spacing and cuddling fixed
=======
void nmod_sparse_mat_mul_vec(mp_ptr y, const nmod_sparse_mat_t M, const mp_ptr x) 
=======
void nmod_sparse_mat_mul_vec(mp_ptr y, const nmod_sparse_mat_t M, mp_srcptr x) 
>>>>>>> Fixed bug in multiply, block Lanczos now works for solving (but not for nullspace)
{
    slong i;
    for (i = 0; i < M->r; ++i) y[i] = nmod_sparse_vec_dot_dense(&M->rows[i], x, M->mod);
>>>>>>> Created sparse vector and matrix utilities for all the fq variants
}
NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_mul_mat(nmod_mat_t Y, const nmod_sparse_mat_t M, const nmod_mat_t X) 
{
    slong i, j;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
    nmod_mat_zero(Y);
    for (i = 0; i < M->r; ++i)
    {
        for (j = 0; j < M->rows[i].nnz; ++j)
        {
            nmod_sparse_entry_struct *e = &M->rows[i].entries[j];
=======
    for (i = 0; i < A->r; ++i)
=======
=======
    nmod_mat_zero(Y);
>>>>>>> Fixed bug in multiply, block Lanczos now works for solving (but not for nullspace)
    for (i = 0; i < M->r; ++i)
>>>>>>> Created sparse vector and matrix utilities for all the fq variants
{
        for (j = 0; j < M->rows[i].nnz; ++j)
{
<<<<<<< HEAD
            nmod_sparse_entry_struct *e = &A->rows[i].entries[j];
>>>>>>> Spacing and cuddling fixed
=======
            nmod_sparse_entry_struct *e = &M->rows[i].entries[j];
>>>>>>> Created sparse vector and matrix utilities for all the fq variants
            _nmod_vec_scalar_addmul_nmod(Y->rows[i], X->rows[e->ind], X->c, e->val, Y->mod);
        }
    }
}

<<<<<<< HEAD
FLINT_DLL
slong nmod_sparse_mat_inv(nmod_sparse_mat_t Ai, const nmod_sparse_mat_t M);
=======
/* Permutations */
/* FLINT_DLL void nmod_sparse_mat_swap_rows(nmod_sparse_mat_t M, slong * perm, slong r, slong s);
FLINT_DLL void nmod_sparse_mat_invert_rows(nmod_sparse_mat_t M, slong * perm);
FLINT_DLL void nmod_sparse_mat_swap_cols(nmod_sparse_mat_t M, slong * perm, slong r, slong s);
FLINT_DLL void nmod_sparse_mat_invert_cols(nmod_sparse_mat_t M, slong * perm);
FLINT_DLL void nmod_sparse_mat_apply_permutation(nmod_sparse_mat_t M, slong * P, slong n);
 */
>>>>>>> Created sparse vector and matrix utilities for all the fq variants

/* Decomposition/reduction */
FLINT_DLL
slong nmod_sparse_mat_lu(slong *P, slong *Q, nmod_sparse_mat_t L, nmod_sparse_mat_t U, const nmod_sparse_mat_t M);
<<<<<<< HEAD

FLINT_DLL
slong nmod_sparse_mat_rref(nmod_sparse_mat_t M);

NMOD_SPARSE_MAT_INLINE
slong nmod_sparse_mat_strong_echelon_form(nmod_sparse_mat_t M)
{
    /* TODO */
    return 0;
}

/* Solve Ax = b */
FLINT_DLL
int nmod_sparse_mat_solve_lanczos(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b, flint_rand_t state);

FLINT_DLL
int nmod_sparse_mat_solve_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b);

FLINT_DLL
int nmod_sparse_mat_solve_lu(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b);

FLINT_DLL
int nmod_sparse_mat_solve_rref(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b);

FLINT_DLL
int nmod_sparse_mat_solve_block_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b, slong block_size, flint_rand_t state);

FLINT_DLL
int nmod_sparse_mat_solve_block_lanczos(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b, slong block_size, flint_rand_t state);

/* Find single nullvector */
FLINT_DLL
int nmod_sparse_mat_nullvector_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, flint_rand_t state);

FLINT_DLL
int nmod_sparse_mat_nullvector_lanczos(mp_ptr x, const nmod_sparse_mat_t M, flint_rand_t state);

FLINT_DLL
int nmod_sparse_mat_nullvector_block_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state); 

FLINT_DLL
int nmod_sparse_mat_nullvector_block_lanczos(mp_ptr x, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state); 

/* Note: this should take in uninitialized matrix X */
FLINT_DLL
slong nmod_sparse_mat_nullspace_lanczos(nmod_mat_t X, const nmod_sparse_mat_t M, flint_rand_t state, slong max_iters);
=======

FLINT_DLL
slong nmod_sparse_mat_rref(nmod_sparse_mat_t M);
>>>>>>> Created sparse vector and matrix utilities for all the fq variants

<<<<<<< HEAD
=======
/* Decomposition/reduction */
>>>>>>> Now with additional utilities, more correct basic functions, and nullspace and inversion functions
FLINT_DLL
slong nmod_sparse_mat_nullspace_block_lanczos(nmod_mat_t X, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state, slong max_iters);

FLINT_DLL
slong nmod_sparse_mat_nullspace_wiedemann(nmod_mat_t X, const nmod_sparse_mat_t M, flint_rand_t state, slong max_iters);

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
/* Solving */
>>>>>>> Now with additional utilities, more correct basic functions, and nullspace and inversion functions
=======
/* Solve Ax=b */
>>>>>>> Added nullvector functions for Lanzcos, everything for basic Wiedemann
=======
/* Solve Ax = b */
>>>>>>> Fixed spacing
FLINT_DLL
<<<<<<< HEAD
slong nmod_sparse_mat_nullspace_block_wiedemann(nmod_mat_t X, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state, slong max_iters);

FLINT_DLL
<<<<<<< HEAD
slong nmod_sparse_mat_nullspace_rref(nmod_mat_t X, const nmod_sparse_mat_t M);
=======
int nmod_sparse_mat_solve_wiedemann(mp_ptr x, const nmod_sparse_mat_t A, const mp_ptr b);

FLINT_DLL
int nmod_sparse_mat_solve_lu(mp_ptr x, const nmod_sparse_mat_t A, const mp_ptr b);
>>>>>>> Added nullvector functions for Lanzcos, everything for basic Wiedemann

FLINT_DLL
slong nmod_sparse_mat_nullspace_lu(nmod_mat_t X, const nmod_sparse_mat_t M);
=======
int nmod_sparse_mat_solve_lanczos(mp_ptr x, const nmod_sparse_mat_t M, const mp_ptr b, flint_rand_t state);

FLINT_DLL
int nmod_sparse_mat_solve_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, const mp_ptr b);

FLINT_DLL
int nmod_sparse_mat_solve_lu(mp_ptr x, const nmod_sparse_mat_t M, const mp_ptr b);

FLINT_DLL
int nmod_sparse_mat_solve_rref(mp_ptr x, const nmod_sparse_mat_t M, const mp_ptr b);

FLINT_DLL
int nmod_sparse_mat_solve_block_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, const mp_ptr b, slong block_size, flint_rand_t state);

FLINT_DLL
int nmod_sparse_mat_solve_block_lanczos(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b, slong block_size, flint_rand_t state);
>>>>>>> Created sparse vector and matrix utilities for all the fq variants

<<<<<<< HEAD
<<<<<<< HEAD
/* Determinant */
FLINT_DLL
mp_limb_t nmod_sparse_mat_det(const nmod_sparse_mat_t M);
=======
=======
/* Find single nullvector */
FLINT_DLL
int nmod_sparse_mat_nullvector_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, flint_rand_t state);

FLINT_DLL
int nmod_sparse_mat_nullvector_lanczos(mp_ptr x, const nmod_sparse_mat_t M, flint_rand_t state);

FLINT_DLL
int nmod_sparse_mat_nullvector_block_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state); 

FLINT_DLL
int nmod_sparse_mat_nullvector_block_lanczos(mp_ptr x, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state); 

>>>>>>> Added nullvector functions for Lanzcos, everything for basic Wiedemann
/* Note: this should take in uninitialized matrix X */
FLINT_DLL
slong nmod_sparse_mat_nullspace_lanczos(nmod_mat_t X, const nmod_sparse_mat_t M, flint_rand_t state, slong max_iters);

FLINT_DLL
slong nmod_sparse_mat_nullspace_block_lanczos(nmod_mat_t X, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state, slong max_iters);

FLINT_DLL
slong nmod_sparse_mat_nullspace_wiedemann(nmod_mat_t X, const nmod_sparse_mat_t M, flint_rand_t state, slong max_iters);

FLINT_DLL
slong nmod_sparse_mat_nullspace_block_wiedemann(nmod_mat_t X, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state, slong max_iters);

FLINT_DLL
slong nmod_sparse_mat_nullspace_rref(nmod_mat_t X, const nmod_sparse_mat_t M);

FLINT_DLL
slong nmod_sparse_mat_nullspace_lu(nmod_mat_t X, const nmod_sparse_mat_t M);

FLINT_DLL
<<<<<<< HEAD
slong nmod_sparse_mat_inv(nmod_sparse_mat_t Ai, const nmod_sparse_mat_t A);
>>>>>>> Now with additional utilities, more correct basic functions, and nullspace and inversion functions
=======
slong nmod_sparse_mat_inv(nmod_sparse_mat_t Ai, const nmod_sparse_mat_t M);
>>>>>>> Created sparse vector and matrix utilities for all the fq variants

/* Nullspace */
/* NMOD_SPARSE_MAT_INLINE
slong nmod_sparse_mat_nullspace(nmod_mat_t X, const nmod_sparse_mat_t M);
 */
#ifdef __cplusplus
}
#endif

#endif

