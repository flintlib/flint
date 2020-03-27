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
#include "nmod_mat.h"
#include "fmpz.h"
#include "thread_support.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Sparse matrix stored in compressed sparse row format */
/* Each entry is a column and a value */
typedef struct 
{
    mp_limb_t val;
    slong col;
} nmod_sparse_mat_entry_struct;

typedef nmod_sparse_mat_entry_struct nmod_sparse_mat_entry_t[1];

/* A sparse matrix is a list of sparse sparse elements */
typedef struct
{
    nmod_sparse_mat_entry_struct *entries;
    slong *row_starts;
    slong *row_nnz;
    slong r;
    slong c;
    slong c_off;
    slong nnz;
    nmod_t mod;
}
nmod_sparse_mat_struct;

typedef nmod_sparse_mat_struct nmod_sparse_mat_t[1];

NMOD_SPARSE_MAT_INLINE
void _nmod_sparse_mat_set_c(nmod_sparse_mat_t mat)
{
    slong i;
    mat->c = 0;
    for(i=0; i<mat->nnz; ++i)
        if(mat->entries[i].col >= mat->c) 
            mat->c = mat->entries[i].col + 1;
}

NMOD_SPARSE_MAT_INLINE
slong nmod_sparse_mat_nrows(const nmod_sparse_mat_t mat)
{
   return mat->r;
}

NMOD_SPARSE_MAT_INLINE
slong nmod_sparse_mat_ncols(const nmod_sparse_mat_t mat)
{
   return mat->c;
}

NMOD_SPARSE_MAT_INLINE
slong nmod_sparse_mat_nnz(const nmod_sparse_mat_t mat)
{
   return mat->nnz;
}

NMOD_SPARSE_MAT_INLINE
void _nmod_sparse_mat_set_mod(nmod_sparse_mat_t mat, mp_limb_t n)
{
    mat->mod.n = n;
    count_leading_zeros(mat->mod.norm, n);
    invert_limb(mat->mod.ninv, n << mat->mod.norm);
}

/* Memory management */
FLINT_DLL void nmod_sparse_mat_init(nmod_sparse_mat_t mat, slong rows, mp_limb_t n);
FLINT_DLL void nmod_sparse_mat_clear(nmod_sparse_mat_t mat);
FLINT_DLL void nmod_sparse_mat_swap(nmod_sparse_mat_t mat1, nmod_sparse_mat_t mat2);

/* One-time instantiation */
FLINT_DLL void nmod_sparse_mat_zero(nmod_sparse_mat_t mat);
FLINT_DLL void nmod_sparse_mat_one(nmod_sparse_mat_t mat);
FLINT_DLL void nmod_sparse_mat_set(nmod_sparse_mat_t mat, const nmod_sparse_mat_t src);
FLINT_DLL void nmod_sparse_mat_set_from_entries(nmod_sparse_mat_t mat, slong * rows, slong * cols, mp_limb_t * vals, slong nnz);

/* Incremental instantiation */
FLINT_DLL void nmod_sparse_mat_append_row(nmod_sparse_mat_t mat, slong row, slong *cols, mp_limb_t *vals, slong nnz);

/* Convert from/to dense matrix */
FLINT_DLL void nmod_sparse_mat_from_dense(nmod_sparse_mat_t mat, const nmod_mat_t src);
FLINT_DLL void nmod_sparse_mat_to_dense(nmod_mat_t mat, const nmod_sparse_mat_t src);

/* Windows and concatenation */
FLINT_DLL void nmod_sparse_mat_window_init(nmod_sparse_mat_t window, const nmod_sparse_mat_t mat, slong r1, slong c1, slong r2, slong c2);
FLINT_DLL void nmod_sparse_mat_window_clear(nmod_sparse_mat_t window);

/* res->r must equal mat1->r + mat2->r, and res->c, mat1->c, and mat2->c must be equal */
FLINT_DLL void nmod_sparse_mat_concat_horizontal(nmod_sparse_mat_t res,
                           const nmod_sparse_mat_t mat1,  const nmod_sparse_mat_t mat2);
/* res->r, mat1->r, and mat2->r must be equal, and res->c must equal mat1->c + mat2->c */
FLINT_DLL void nmod_sparse_mat_concat_vertical(nmod_sparse_mat_t res,
                           const nmod_sparse_mat_t mat1,  const nmod_sparse_mat_t mat2);

/* Random matrix generation */
FLINT_DLL void nmod_sparse_mat_randtest(nmod_sparse_mat_t mat, flint_rand_t state);
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

FLINT_DLL int nmod_sparse_mat_equal(const nmod_sparse_mat_t mat1, const nmod_sparse_mat_t mat2);

NMOD_SPARSE_MAT_INLINE
int nmod_sparse_mat_is_zero(const nmod_sparse_mat_t mat) {
    return mat->nnz == 0;
}

NMOD_SPARSE_MAT_INLINE
int
nmod_sparse_mat_is_zero_row(const nmod_sparse_mat_t mat, slong i)
{
    return mat->row_nnz[i]==0;
}

NMOD_SPARSE_MAT_INLINE
int nmod_sparse_mat_is_square(const nmod_sparse_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Must have A->r == B->c and A->c == B->r */
FLINT_DLL void nmod_sparse_mat_transpose(nmod_sparse_mat_t B, const nmod_sparse_mat_t A);

/* Addition and subtraction */
/* Arguments must be distinct */
FLINT_DLL void nmod_sparse_mat_add(nmod_sparse_mat_t C, const nmod_sparse_mat_t A, const nmod_sparse_mat_t B);
FLINT_DLL void nmod_sparse_mat_sub(nmod_sparse_mat_t C, const nmod_sparse_mat_t A, const nmod_sparse_mat_t B);

FLINT_DLL void nmod_sparse_mat_neg(nmod_sparse_mat_t B, const nmod_sparse_mat_t A);

/* Matrix-scalar arithmetic */

FLINT_DLL void nmod_sparse_mat_scalar_mul(nmod_sparse_mat_t B, const nmod_sparse_mat_t A, mp_limb_t c);

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_scalar_mul_fmpz(nmod_sparse_mat_t res, const nmod_sparse_mat_t M, const fmpz_t c)
{
    fmpz_t d;
    fmpz_init(d);
    fmpz_mod_ui(d, c, res->mod.n);
    nmod_sparse_mat_scalar_mul(res, M, fmpz_get_ui(d));
    fmpz_clear(d);
}

/* Matrix-vector and matrix-matrix multipliciation */
FLINT_DLL void nmod_sparse_mat_mul_vec(mp_ptr y, const nmod_sparse_mat_t A, const mp_ptr x);
FLINT_DLL void nmod_sparse_mat_mul_mat(nmod_mat_t Y, const nmod_sparse_mat_t A, const nmod_mat_t X);

/* Permutations */
/* FLINT_DLL void nmod_sparse_mat_swap_rows(nmod_sparse_mat_t mat, slong * perm, slong r, slong s);
FLINT_DLL void nmod_sparse_mat_invert_rows(nmod_sparse_mat_t mat, slong * perm);
FLINT_DLL void nmod_sparse_mat_swap_cols(nmod_sparse_mat_t mat, slong * perm, slong r, slong s);
FLINT_DLL void nmod_sparse_mat_invert_cols(nmod_sparse_mat_t mat, slong * perm);
FLINT_DLL void nmod_sparse_mat_apply_permutation(nmod_sparse_mat_t A, slong * P, slong n);
 */
/* Nonsingular solving */
int nmod_sparse_mat_solve_lanczos(mp_ptr x, const nmod_sparse_mat_t A, const mp_ptr b, flint_rand_t state);

/* FLINT_DLL int nmod_sparse_mat_solve(nmod_mat_t X, const nmod_sparse_mat_t A, const nmod_mat_t B);
FLINT_DLL int nmod_sparse_mat_solve_vec(mp_ptr x, const nmod_sparse_mat_t A, mp_srcptr b);
 */
/* Nullspace */

/* FLINT_DLL slong nmod_mat_nullspace(nmod_mat_t X, const nmod_sparse_mat_t A);
 */
/*
   Suggested initial modulus size for multimodular algorithms. This should
   be chosen so that we get the most number of bits per cycle
   in matrix multiplication. On x86-64 it appears to be optimal to use
   moduli giving nlimbs = 2. This should hold both in the classical
   range and in Strassen blocks.
 */
#define NMOD_SPARSE_MAT_OPTIMAL_MODULUS_BITS (FLINT_BITS-5)

#ifdef __cplusplus
}
#endif

#endif

