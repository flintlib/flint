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
#ifndef FMPZ_SPARSE_MAT_H
#define FMPZ_SPARSE_MAT_H

#ifdef FMPZ_SPARSE_MAT_INLINES_C
#define FMPZ_SPARSE_MAT_INLINE FLINT_DLL
#else
#define FMPZ_SPARSE_MAT_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "hashmap.h"
#include "heap.h"
#include "nmod_sparse_mat.h"
#include "fmpz_sparse_vec.h"
#include "fmpz_mat.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* A sparse matrix is a list of sparse row vectors */
typedef struct
{
    fmpz_sparse_vec_struct *rows;
    slong r;
    slong c;
    slong c_off;
}
fmpz_sparse_mat_struct;

typedef fmpz_sparse_mat_struct fmpz_sparse_mat_t[1];

#define LT(M, r) M->rows[r].entries[0]

/* Can package matrix with its column support, to enable various operations */
typedef struct
{
    fmpz_sparse_mat_struct *M;
    hashmap_struct *cols;
}
fmpz_sparse_mat_with_transpose_struct;

typedef fmpz_sparse_mat_with_transpose_struct fmpz_sparse_mat_with_transpose_t[1];

/* Memory management */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_init (fmpz_sparse_mat_t M, slong rows, slong cols) 
{
    M->rows = flint_calloc(rows, sizeof(*M->rows));
    M->r = rows;
    M->c = cols;
    M->c_off = 0;
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_clear (fmpz_sparse_mat_t M) 
{
    slong i;
    for (i = 0; i < M->r; ++i) fmpz_sparse_vec_clear(&M->rows[i]);
    flint_free(M->rows);
    memset(M, 0, sizeof(*M));
}
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_swap (fmpz_sparse_mat_t M1, fmpz_sparse_mat_t M2) 
{
    fmpz_sparse_mat_t tmp;
    *tmp = *M1; *M1 = *M2; *M2 = *tmp;
}

FLINT_DLL
slong fmpz_sparse_mat_max_bits(const fmpz_sparse_mat_t v);

/* One-time instantiation */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_zero (fmpz_sparse_mat_t M) 
{
    slong i;
    for (i = 0; i < M->r; ++i) fmpz_sparse_vec_zero(&M->rows[i]);
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_one (fmpz_sparse_mat_t M) 
{
    slong i;
    for (i = 0; i < M->r; ++i) fmpz_sparse_vec_one(&M->rows[i], i);
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_set (fmpz_sparse_mat_t N, const fmpz_sparse_mat_t M) 
{
    slong i, rmax = FLINT_MIN(M->r, M->r);
    if (M==N) return;
    for (i = 0; i < rmax; ++i) fmpz_sparse_vec_set(&N->rows[i], &M->rows[i], M->c_off);
}

FLINT_DLL
void fmpz_sparse_mat_from_entries(fmpz_sparse_mat_t M, slong * rows, slong * cols, fmpz * vals, slong nnz);

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_append_col (fmpz_sparse_mat_t M, const fmpz *v) 
{
    slong i;
    for (i = 0; i < M->r; ++i) fmpz_sparse_vec_set_entry(&M->rows[i], M->c, &v[i]);
    M->c += 1;
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_append_row (fmpz_sparse_mat_t M, const fmpz_sparse_vec_t v) 
{
    M->rows = realloc(M->rows, (M->r+1)*sizeof(*M->rows));
    memset(M->rows + M->r, 0, sizeof(*M->rows));
    fmpz_sparse_vec_set(&M->rows[M->r], v, 0);
    M->r += 1;
}

/* Convert from/to dense matrix */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_from_dense (fmpz_sparse_mat_t M, const fmpz_mat_t dM)
{
    slong i, rmax = FLINT_MIN(M->r, dM->r);
    if (M->c == 0 || dM->c == 0) {fmpz_sparse_mat_zero(M); return;}
    for (i = 0; i < rmax; ++i) fmpz_sparse_vec_from_dense(&M->rows[i], dM->rows[i], dM->c);
}
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_to_dense (fmpz_mat_t dM, const fmpz_sparse_mat_t M)
{
    slong i, rmax = FLINT_MIN(M->r, dM->r);
    if (M->c == 0 || dM->c == 0) {fmpz_mat_zero(dM); return;}
    for (i = 0; i < rmax; ++i) 
        fmpz_sparse_vec_to_dense(dM->rows[i], &M->rows[i], dM->c);
}

/* Convert to/from modular matrix */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_get_nmod_sparse_mat(nmod_sparse_mat_t Amod, const fmpz_sparse_mat_t A)
{
    slong i;
    for (i = 0; i < A->r; i++)
        fmpz_sparse_vec_get_nmod_sparse_vec(&Amod->rows[i], &A->rows[i], Amod->mod);
}

void fmpz_sparse_mat_multi_mod_ui_precomp(nmod_sparse_mat_struct * residues, slong nres, const fmpz_sparse_mat_t M, 
                                          const fmpz_comb_t comb, fmpz_comb_temp_t temp);

void fmpz_sparse_mat_multi_mod_ui(nmod_sparse_mat_struct * residues, slong nres, const fmpz_sparse_mat_t M);

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_set_nmod_sparse_mat_unsigned(fmpz_sparse_mat_t A, const nmod_sparse_mat_t Amod)
{
    slong i;
    for (i = 0; i < A->r; i++)
        fmpz_sparse_vec_set_nmod_sparse_vec_unsigned(&A->rows[i], &Amod->rows[i]);
}


FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_set_nmod_sparse_mat(fmpz_sparse_mat_t A, const nmod_sparse_mat_t Amod)
{
    slong i;
    for (i = 0; i < A->r; i++)
        fmpz_sparse_vec_set_nmod_sparse_vec(&A->rows[i], &Amod->rows[i], Amod->mod);
}

FLINT_DLL
void fmpz_sparse_mat_CRT_ui(fmpz_sparse_mat_t res, const fmpz_sparse_mat_t mat1,
                        const fmpz_t m1, const nmod_sparse_mat_t mat2, int sign);


FLINT_DLL
void fmpz_sparse_mat_multi_CRT_ui_precomp(fmpz_sparse_mat_t M, nmod_sparse_mat_struct * residues, slong nres,
                                          const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign);

FLINT_DLL
void fmpz_sparse_mat_multi_CRT_ui(fmpz_sparse_mat_t mat, nmod_sparse_mat_struct * residues, slong nres, int sign);

/* Windows, concatenation, and splitting */
FLINT_DLL
void fmpz_sparse_mat_window_init (fmpz_sparse_mat_t W, const fmpz_sparse_mat_t M, slong r1, slong c1, slong r2, slong c2);

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_window_clear (fmpz_sparse_mat_t W) 
{
    flint_free(W->rows);
    memset(W, 0, sizeof(*W));
}


/* Combine M1 and M2 into block matrix B = [M1 M2] */
/* B->r must equal M1->r and M2->r */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_concat_horizontal(fmpz_sparse_mat_t B,
                                    const fmpz_sparse_mat_t M1,  const fmpz_sparse_mat_t M2) 
{
    slong i;
    B->c = M1->c + M2->c;
    for (i = 0; i < B->r; ++i)
        fmpz_sparse_vec_concat(&B->rows[i], &M1->rows[i], &M2->rows[i], M1->c);
}
/* Combine M1 and M2 into block matrix B = [M1^t M1^t]^t */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_concat_vertical (fmpz_sparse_mat_t B, const fmpz_sparse_mat_t M1,  const fmpz_sparse_mat_t M2) 
{
    slong i;
    B->c = FLINT_MAX(M1->c, M2->c);
    for (i = 0; i < M1->r; ++i)
        fmpz_sparse_vec_set(&B->rows[i], &M1->rows[i], M1->c_off);
    for (i = M1->r; i < B->r; ++i)
        fmpz_sparse_vec_set(&B->rows[i], &M2->rows[i-M1->r], M2->c_off);
}

/* Split block matrix B = [M1 M2] into submatrices M1 and M2 */
/* M1->r and M2->r must equal B->r */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_split_horizontal (fmpz_sparse_mat_t M1, fmpz_sparse_mat_t M2, const fmpz_sparse_mat_t B, slong c)
{
    slong i;
    for (i = 0; i < B->r; ++i) fmpz_sparse_vec_split(&M1->rows[i], &M2->rows[i], &B->rows[i], c);
}

/* Split block matix B = [M1^t M1^t]^t into submatrices M1 and M2 */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_split_vertical (fmpz_sparse_mat_t M1, fmpz_sparse_mat_t M2, const fmpz_sparse_mat_t B, slong r)
{
    slong i;
    r = FLINT_MIN(r, B->r);
    for (i = 0; i < r; ++i) fmpz_sparse_vec_set(&M1->rows[i], &B->rows[i], B->c_off);
    for (i = r; i < B->r; ++i) fmpz_sparse_vec_set(&M2->rows[i-r], &B->rows[i], B->c_off);
}

/* Matrix permutation */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_permute_cols(fmpz_sparse_mat_t M, slong *Q) 
{
    slong i;
    for (i = 0; i < M->r; ++i) 
    {
        if (!M->rows[i].nnz) continue;
        fmpz_sparse_vec_permute_inds(&M->rows[i], Q);
        qsort(M->rows[i].entries, M->rows[i].nnz, sizeof(*M->rows[i].entries), fmpz_sparse_entry_cmp);
    }
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_permute_rows(fmpz_sparse_mat_t M, slong *P) 
{
    slong i;
    fmpz_sparse_vec_struct *prows;
    prows = flint_calloc(M->r, sizeof(*prows));
    for (i = 0; i < M->r; ++i) prows[P[i]] = M->rows[i];
    memcpy(M->rows, prows, M->r*sizeof(*M->rows));
    flint_free(prows);
}

/* Random matrix generation */
FLINT_DLL void fmpz_sparse_mat_randtest (fmpz_sparse_mat_t M, flint_rand_t state, slong min_nnz, slong max_nnz, flint_bitcnt_t bits);
FLINT_DLL void fmpz_sparse_mat_randtest_unsigned (fmpz_sparse_mat_t M, flint_rand_t state, slong min_nnz, slong max_nnz, flint_bitcnt_t bits);
/*
FLINT_DLL void fmpz_sparse_mat_randfull (fmpz_sparse_mat_t M, flint_rand_t state, fmpz_ctx_t ctx);
FLINT_DLL int fmpz_sparse_mat_randpermdiag(fmpz_sparse_mat_t M, flint_rand_t state,
                 const fmpz *diag, slong n);
FLINT_DLL void fmpz_sparse_mat_randrank (fmpz_sparse_mat_t, flint_rand_t state, slong rank, fmpz_ctx_t ctx);
FLINT_DLL void fmpz_sparse_mat_randops (fmpz_sparse_mat_t M, slong count, flint_rand_t state, fmpz_ctx_t ctx);
FLINT_DLL void fmpz_sparse_mat_randtril (fmpz_sparse_mat_t M, flint_rand_t state, int unit, fmpz_ctx_t ctx);
FLINT_DLL void fmpz_sparse_mat_randtriu (fmpz_sparse_mat_t M, flint_rand_t state, int unit, fmpz_ctx_t ctx);
 */

FLINT_DLL void fmpz_sparse_mat_print_pretty (const fmpz_sparse_mat_t M);

FMPZ_SPARSE_MAT_INLINE
int fmpz_sparse_mat_equal (const fmpz_sparse_mat_t M1, const fmpz_sparse_mat_t M2) 
{
    slong i;
    if (M1->r != M2->r) return 0;
    for (i = 0; i < M1->r; ++i)
        if (!fmpz_sparse_vec_equal(&M1->rows[i], &M2->rows[i], M1->c_off-M2->c_off)) return 0;
    return 1;
}

FMPZ_SPARSE_MAT_INLINE
int fmpz_sparse_mat_is_zero (const fmpz_sparse_mat_t M) 
{
    slong i;
    for (i = 0; i < M->r; ++i) 
        if (!fmpz_sparse_vec_is_zero(&M->rows[i])) return 0;
    return 1;
}

/* Must have M->r == N->c and M->c == N->r */
FLINT_DLL 
void fmpz_sparse_mat_transpose (fmpz_sparse_mat_t N, const fmpz_sparse_mat_t M);

/* Arithmetic */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_neg (fmpz_sparse_mat_t N, const fmpz_sparse_mat_t M) 
{
    slong i;
    for (i = 0; i < N->r; ++i) fmpz_sparse_vec_neg(&N->rows[i], &M->rows[i]);
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_scalar_mul_fmpz(fmpz_sparse_mat_t N, const fmpz_sparse_mat_t M, const fmpz_t c) 
{
    slong i;
    for (i = 0; i < N->r; ++i) fmpz_sparse_vec_scalar_mul_fmpz(&N->rows[i], &M->rows[i], c);    
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_mul_diag_fmpz(fmpz_sparse_mat_t N, const fmpz_sparse_mat_t M, fmpz * D) 
{
    slong i;
    for (i = 0; i < N->r; ++i) fmpz_sparse_vec_scalar_mul_fmpz(&N->rows[i], &M->rows[i], &D[i]);    
}


FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_scalar_divexact_fmpz(fmpz_sparse_mat_t N, const fmpz_sparse_mat_t M, const fmpz_t c) 
{
    slong i;
    for (i = 0; i < N->r; ++i) fmpz_sparse_vec_scalar_divexact_fmpz(&N->rows[i], &M->rows[i], c);    
}


FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_scalar_mod_fmpz(fmpz_sparse_mat_t N, const fmpz_sparse_mat_t M, const fmpz_t mod) 
{
    slong i;
    for (i = 0; i < N->r; ++i) fmpz_sparse_vec_scalar_mod_fmpz(&N->rows[i], &M->rows[i], mod);    
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_scalar_mods_fmpz(fmpz_sparse_mat_t N, const fmpz_sparse_mat_t M, const fmpz_t mod) 
{
    slong i;
    for (i = 0; i < N->r; ++i) fmpz_sparse_vec_scalar_mods_fmpz(&N->rows[i], &M->rows[i], mod);    
}


FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_add (fmpz_sparse_mat_t O, const fmpz_sparse_mat_t M, const fmpz_sparse_mat_t N) 
{
    slong i;
    for (i = 0; i < O->r; ++i) fmpz_sparse_vec_add(&O->rows[i], &M->rows[i], &N->rows[i]);
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_sub (fmpz_sparse_mat_t O, const fmpz_sparse_mat_t M, const fmpz_sparse_mat_t N) 
{
    slong i;
    for (i = 0; i < O->r; ++i) fmpz_sparse_vec_sub(&O->rows[i], &M->rows[i], &N->rows[i]);
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_scalar_addmul_fmpz(fmpz_sparse_mat_t O, const fmpz_sparse_mat_t M, const fmpz_sparse_mat_t N, const fmpz_t c) 
{
    slong i;
    for (i = 0; i < O->r; ++i) fmpz_sparse_vec_scalar_addmul_fmpz(&O->rows[i], &M->rows[i], &N->rows[i], c);
}

FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_scalar_submul_fmpz(fmpz_sparse_mat_t O, const fmpz_sparse_mat_t M, const fmpz_sparse_mat_t N, const fmpz_t c) 
{
    slong i;
    for (i = 0; i < O->r; ++i) fmpz_sparse_vec_scalar_submul_fmpz(&O->rows[i], &M->rows[i], &N->rows[i], c);
}

/* Matrix-vector and matrix-matrix multipliciation */
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_mul_vec (fmpz *y, const fmpz_sparse_mat_t M, const fmpz *x) 
{
    slong i;
    if (M->c == 0) _fmpz_vec_zero(y, M->r);
    else for (i = 0; i < M->r; ++i) fmpz_sparse_vec_dot_dense(&y[i], &M->rows[i], x);
}
FMPZ_SPARSE_MAT_INLINE
void fmpz_sparse_mat_mul_mat (fmpz_mat_t Y, const fmpz_sparse_mat_t M, const fmpz_mat_t X) 
{
    slong i, j;
    fmpz_mat_zero (Y);
    for (i = 0; i < M->r; ++i)
    {
        for (j = 0; j < M->rows[i].nnz; ++j)
        {
            fmpz_sparse_entry_struct *e = &M->rows[i].entries[j];
            _fmpz_vec_scalar_addmul_fmpz(Y->rows[i], X->rows[e->ind], X->c, e->val);
        }
    }
}

/* Memory management for matrix with transpose */
FLINT_DLL
void _fmpz_sparse_mat_with_transpose_init(fmpz_sparse_mat_with_transpose_t MT, fmpz_sparse_mat_t M);

FMPZ_SPARSE_MAT_INLINE
void _fmpz_sparse_mat_with_transpose_clear(fmpz_sparse_mat_with_transpose_t MT)
{
    slong c;
    for (c = 0; c < MT->M->c; ++c)
        hashmap_clear(&MT->cols[c]);
    flint_free(MT->cols);
    memset(MT, 0, sizeof(*MT));
}

FMPZ_SPARSE_MAT_INLINE
void _fmpz_sparse_mat_with_transpose_print_pretty(fmpz_sparse_mat_with_transpose_t MT)
{
    slong i, j;
    fmpz_sparse_mat_print_pretty(MT->M);
    flint_printf("Transpose: \nSupport: ");
    for (i = 0; i < MT->M->c; ++i)
    {
        flint_printf("%wd ", MT->cols[i].num);
    }
    flint_printf("\n");
    for (i = 0; i < MT->M->c; ++i)
    {
        flint_printf("%wd: [", i); fflush(stdout);
        for (j = 0; j < MT->cols[i].num; ++j)
        {
            if (j > 0) flint_printf(" ");
            flint_printf("%wd: ", MT->cols[i].keys[j]);  fflush(stdout);
            fmpz_print(*((fmpz_t *) (MT->cols[i].vals[j])));
        }
        flint_printf("]\n");
    }
}

FLINT_DLL
void _fmpz_sparse_mat_with_transpose_fix_support(fmpz_sparse_mat_with_transpose_t MT, slong r, slong *osupp, slong onnz);

#define MT_FIX(MT, r, ...) \
{\
    slong *supp, nnz; \
    nnz = _fmpz_sparse_vec_support(&supp, &MT->M->rows[r]); \
    __VA_ARGS__; \
    _fmpz_sparse_mat_with_transpose_fix_support(MT, r, supp, nnz);\
    flint_free(supp);\
}

FMPZ_SPARSE_MAT_INLINE
int _fmpz_sparse_mat_with_transpose_gauss_elim_col(fmpz_sparse_mat_with_transpose_t MT, slong pr, slong r, slong col)
{
    MT_FIX(MT, r, fmpz_sparse_vec_gauss_elim_col(&MT->M->rows[r], &MT->M->rows[pr], col));
    return fmpz_sparse_vec_at(&MT->M->rows[r], col) == NULL;
}

FMPZ_SPARSE_MAT_INLINE
int _fmpz_sparse_mat_with_transpose_gauss_elim(fmpz_sparse_mat_with_transpose_t MT, slong pr, slong r)
{
    MT_FIX(MT, r, fmpz_sparse_vec_gauss_elim(&MT->M->rows[r], &MT->M->rows[pr]));
    return fmpz_sparse_vec_at(&MT->M->rows[r], MT->M->rows[pr].entries[0].ind) == NULL;
}

FMPZ_SPARSE_MAT_INLINE
int _fmpz_sparse_mat_with_transpose_gauss_elim_mod(fmpz_sparse_mat_with_transpose_t MT, slong pr, slong r, const fmpz_t mod)
{
    MT_FIX(MT, r, 
    fmpz_sparse_vec_gauss_elim(&MT->M->rows[r], &MT->M->rows[pr]);
    fmpz_sparse_vec_scalar_mod_fmpz(&MT->M->rows[r], &MT->M->rows[r], mod));
    return fmpz_sparse_vec_at(&MT->M->rows[r], MT->M->rows[pr].entries[0].ind) == NULL;
}

FMPZ_SPARSE_MAT_INLINE
int _fmpz_sparse_mat_with_transpose_gauss_elim_mods(fmpz_sparse_mat_with_transpose_t MT, slong pr, slong r, const fmpz_t mod)
{
    MT_FIX(MT, r, 
    fmpz_sparse_vec_gauss_elim(&MT->M->rows[r], &MT->M->rows[pr]);
    fmpz_sparse_vec_scalar_mods_fmpz(&MT->M->rows[r], &MT->M->rows[r], mod));
    return fmpz_sparse_vec_at(&MT->M->rows[r], MT->M->rows[pr].entries[0].ind) == NULL;
}

FMPZ_SPARSE_MAT_INLINE
int _fmpz_sparse_mat_with_transpose_gauss_elim_ext(fmpz_sparse_mat_with_transpose_t MT, slong pr, slong r)
{
    /* If leading entries do not match, or leading entry of pr divides that of r, just a normal elimination */
    fmpz_sparse_entry_struct *pe = &MT->M->rows[pr].entries[0], *e = &MT->M->rows[r].entries[0]; 
    if (pe->ind != e->ind)
        return _fmpz_sparse_mat_with_transpose_gauss_elim(MT, pr, r);
    
    MT_FIX(MT, pr, MT_FIX(MT, r, 
    fmpz_sparse_vec_gauss_elim_ext(&MT->M->rows[r], &MT->M->rows[pr])
    ));
    return 1;
}

FMPZ_SPARSE_MAT_INLINE
int _fmpz_sparse_mat_with_transpose_gauss_elim_ext_mod(fmpz_sparse_mat_with_transpose_t MT, slong pr, slong r, const fmpz_t mod)
{
    /* If leading entries do not match, or leading entry of pr divides that of r, just a normal elimination */
    fmpz_sparse_entry_struct *pe = &MT->M->rows[pr].entries[0], *e = &MT->M->rows[r].entries[0]; 
    if (pe->ind != e->ind)
        return _fmpz_sparse_mat_with_transpose_gauss_elim_mod(MT, pr, r, mod);
    
    MT_FIX(MT, pr, MT_FIX(MT, r, 
    fmpz_sparse_vec_gauss_elim_ext(&MT->M->rows[r], &MT->M->rows[pr]);
    fmpz_sparse_vec_scalar_mod_fmpz(&MT->M->rows[pr], &MT->M->rows[pr], mod);
    fmpz_sparse_vec_scalar_mod_fmpz(&MT->M->rows[r], &MT->M->rows[r], mod);
    ));
    return 1;
}

FMPZ_SPARSE_MAT_INLINE
int _fmpz_sparse_mat_with_transpose_gauss_elim_ext_mods(fmpz_sparse_mat_with_transpose_t MT, slong pr, slong r, const fmpz_t mod)
{
    /* If leading entries do not match, or leading entry of pr divides that of r, just a normal elimination */
    fmpz_sparse_entry_struct *pe = &MT->M->rows[pr].entries[0], *e = &MT->M->rows[r].entries[0]; 
    if (pe->ind != e->ind)
        return _fmpz_sparse_mat_with_transpose_gauss_elim_mods(MT, pr, r, mod);
    
    MT_FIX(MT, pr, MT_FIX(MT, r,
    fmpz_sparse_vec_gauss_elim_ext(&MT->M->rows[r], &MT->M->rows[pr]);
    fmpz_sparse_vec_scalar_mods_fmpz(&MT->M->rows[pr], &MT->M->rows[pr], mod);
    fmpz_sparse_vec_scalar_mods_fmpz(&MT->M->rows[r], &MT->M->rows[r], mod);
    ));
    return 1;
}

/* Utility computations */
FLINT_DLL
void fmpz_sparse_mat_content(fmpz_t mat_gcd, const fmpz_sparse_mat_t M);

FLINT_DLL
void fmpz_sparse_mat_gram(fmpz_mat_t B, const fmpz_sparse_mat_t A);

/* Solving */
FLINT_DLL
void fmpz_sparse_mat_solve_bound(fmpz_t N, fmpz_t D, const fmpz_sparse_mat_t A, const fmpz_mat_t B);

FLINT_DLL
int fmpz_sparse_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod, const fmpz_sparse_mat_t A, const fmpz_mat_t B);

FMPZ_SPARSE_MAT_INLINE
int fmpz_sparse_mat_solve_vec_dixon(fmpz * x, fmpz_t mod, const fmpz_sparse_mat_t A, fmpz *b)
{
    int ret;
    slong i;
    fmpz_mat_t X, B;
    fmpz_mat_init(X, A->c, 1);
    fmpz_mat_init(B, A->r, 1);
    for (i = 0; i < A->r; ++i) fmpz_set(fmpz_mat_entry(B, i, 0), &b[i]);
    ret = fmpz_sparse_mat_solve_dixon(X, mod, A, B);
    for (i = 0; i < A->c; ++i) fmpz_set(&x[i], fmpz_mat_entry(X, i, 0));
    fmpz_mat_clear(X);
    fmpz_mat_clear(B);
    return ret;
}

FLINT_DLL
int fmpz_sparse_mat_solve_dixon_den(fmpz_mat_t X, fmpz_t den, const fmpz_sparse_mat_t A, const fmpz_mat_t B);

FMPZ_SPARSE_MAT_INLINE
int fmpz_sparse_mat_solve_vec_dixon_den(fmpz * x, fmpz_t den, const fmpz_sparse_mat_t A, fmpz *b)
{
    int ret;
    slong i;
    fmpz_mat_t X, B;
    fmpz_mat_init(X, A->c, 1);
    fmpz_mat_init(B, A->r, 1);
    for (i = 0; i < A->r; ++i) fmpz_set(fmpz_mat_entry(B, i, 0), &b[i]);
    ret = fmpz_sparse_mat_solve_dixon_den(X, den, A, B);
    for (i = 0; i < A->c; ++i) fmpz_set(&x[i], fmpz_mat_entry(X, i, 0));
    fmpz_mat_clear(X);
    fmpz_mat_clear(B);
    return ret;
}

/* Determinant computation */
FLINT_DLL
void fmpz_sparse_mat_det_bound(fmpz_t bound, const fmpz_sparse_mat_t A);
FLINT_DLL
void fmpz_sparse_mat_det_cofactor(fmpz_t det, const fmpz_sparse_mat_t M);
FLINT_DLL
void fmpz_sparse_mat_det_bareiss(fmpz_t det, const fmpz_sparse_mat_t M);
FLINT_DLL
void fmpz_sparse_mat_det_divisor(fmpz_t d, const fmpz_sparse_mat_t M);
FLINT_DLL
void fmpz_sparse_mat_det_modular_accelerated(fmpz_t det, const fmpz_sparse_mat_t A, int proved);
FLINT_DLL
void fmpz_sparse_mat_det_modular_given_divisor(fmpz_t det, const fmpz_sparse_mat_t A, const fmpz_t d, int proved);
FLINT_DLL
void fmpz_sparse_mat_det_modular(fmpz_t det, const fmpz_sparse_mat_t A, int proved);
FLINT_DLL
void fmpz_sparse_mat_det(fmpz_t det, const fmpz_sparse_mat_t A);

/* Fraction-free LU factorization */
FLINT_DLL
slong fmpz_sparse_mat_fflu(fmpz *D, slong *P, slong *Q, fmpz_sparse_mat_t L, fmpz_sparse_mat_t U, 
                            const fmpz_sparse_mat_t M);

/* Hermite normal form */
FLINT_DLL
int fmpz_sparse_mat_is_in_hnf(const fmpz_sparse_mat_t A);

FLINT_DLL
slong fmpz_sparse_mat_hnf_classical(fmpz_sparse_mat_t M);

FLINT_DLL
slong fmpz_sparse_mat_hnf_xgcd(fmpz_sparse_mat_t M);

FLINT_DLL
slong fmpz_sparse_mat_hnf_minors(fmpz_sparse_mat_t M);

FLINT_DLL
slong fmpz_sparse_mat_hnf_modular(fmpz_sparse_mat_t M, const fmpz_t det);

FLINT_DLL
slong fmpz_sparse_mat_hnf_modular_eldiv(fmpz_sparse_mat_t M, const fmpz_t n);

/* Modular forms */
slong fmpz_sparse_mat_howell_form_mod(fmpz_sparse_mat_t M, const fmpz_t mod);

FLINT_DLL
slong fmpz_sparse_mat_strong_echelon_form_mod(fmpz_sparse_mat_t M, const fmpz_t mod);

#ifdef __cplusplus
}
#endif

#endif

