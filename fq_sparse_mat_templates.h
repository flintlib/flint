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

#ifdef __cplusplus
 extern "C" {
#endif

/* A sparse matrix is a list of sparse vectors */
typedef struct
{
    TEMPLATE(T, sparse_vec_struct) *rows;
    slong r;
    slong c;
    slong c_off;
}
TEMPLATE(T, sparse_mat_struct);

typedef TEMPLATE(T, sparse_mat_struct) TEMPLATE(T, sparse_mat_t)[1];

/* Memory management */
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_init) (TEMPLATE(T, sparse_mat_t) M, slong rows, slong cols, const TEMPLATE(T, ctx_t) ctx) 
{
    M->rows = flint_calloc(rows, sizeof(*M->rows));
    M->r = rows;
    M->c = cols;
    M->c_off = 0;
}
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_clear) (TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < M->r; ++i) TEMPLATE(T, sparse_vec_clear)(&M->rows[i], ctx);
    flint_free(M->rows);
    memset(M, 0, sizeof(*M));
}
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_swap) (TEMPLATE(T, sparse_mat_t) M1, TEMPLATE(T, sparse_mat_t) M2, const TEMPLATE(T, ctx_t) ctx) 
{
    TEMPLATE(T, sparse_mat_t) tmp;
    *tmp = *M1; *M1 = *M2; *M2 = *tmp;
}

/* One-time instantiation */
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_zero) (TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < M->r; ++i) TEMPLATE(T, sparse_vec_zero)(&M->rows[i], ctx);
}

FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_one) (TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < M->r; ++i) TEMPLATE(T, sparse_vec_one)(&M->rows[i], i, ctx);
}

FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_set) (TEMPLATE(T, sparse_mat_t) N, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i, rmax = FLINT_MIN(M->r, M->r);
    if(M==N) return;
    for(i=0; i<rmax; ++i) TEMPLATE(T, sparse_vec_set)(&N->rows[i], &M->rows[i], M->c_off, ctx);
}

FLINT_DLL
void TEMPLATE(T, sparse_mat_from_entries)(TEMPLATE(T, sparse_mat_t) M, slong * rows, slong * cols, TEMPLATE(T, struct) * vals, slong nnz, const TEMPLATE(T, ctx_t) ctx);

FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_append_col) (TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *v, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for(i=0; i<M->r; ++i) TEMPLATE(T, sparse_vec_set_entry)(&M->rows[i], M->c, &v[i], ctx);
    M->c += 1;
}

FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_append_row) (TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, sparse_vec_t) v, const TEMPLATE(T, ctx_t) ctx) 
{
    M->rows = realloc(M->rows, (M->r+1)*sizeof(*M->rows));
    memset(M->rows + M->r, 0, sizeof(*M->rows));
    TEMPLATE(T, sparse_vec_set)(&M->rows[M->r], v, 0, ctx);
    M->r += 1;
}

/* Convert from/to dense matrix */
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_from_dense) (TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, mat_t) dM, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, rmax = FLINT_MIN(M->r, dM->r);
    for (i = 0; i < rmax; ++i) TEMPLATE(T, sparse_vec_from_dense)(&M->rows[i], dM->rows[i], dM->c, ctx);
}
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_to_dense) (TEMPLATE(T, mat_t) dM, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, rmax = FLINT_MIN(M->r, dM->r);
    for (i = 0; i < rmax; ++i) TEMPLATE(T, sparse_vec_to_dense)(dM->rows[i], &M->rows[i], dM->c, ctx);
}

/* Windows, concatenation, and splitting */
FLINT_DLL
void TEMPLATE(T, sparse_mat_window_init) (TEMPLATE(T, sparse_mat_t) W, const TEMPLATE(T, sparse_mat_t) M, slong r1, slong c1, slong r2, slong c2, const TEMPLATE(T, ctx_t) ctx);

FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_window_clear) (TEMPLATE(T, sparse_mat_t) W, TEMPLATE(T, ctx_t) ctx) 
{
    flint_free(W->rows);
    memset(W, 0, sizeof(*W));
}


/* Combine M1 and M2 into block matrix B = [M1 M2] */
/* B->r must equal M1->r and M2->r */
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_concat_horizontal)(TEMPLATE(T, sparse_mat_t) B,
                                    const TEMPLATE(T, sparse_mat_t) M1,  const TEMPLATE(T, sparse_mat_t) M2, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    B->c = M1->c + M2->c;
    for (i = 0; i < B->r; ++i)
        TEMPLATE(T, sparse_vec_concat)(&B->rows[i], &M1->rows[i], &M2->rows[i], M1->c, ctx);
}
/* Combine M1 and M2 into block matrix B = [M1^t M1^t]^t */
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_concat_vertical) (TEMPLATE(T, sparse_mat_t) B, const TEMPLATE(T, sparse_mat_t) M1,  const TEMPLATE(T, sparse_mat_t) M2, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    B->c = FLINT_MAX(M1->c, M2->c);
    for (i = 0; i < M1->r; ++i)
        TEMPLATE(T, sparse_vec_set)(&B->rows[i], &M1->rows[i], M1->c_off, ctx);
    for (i = M1->r; i < B->r; ++i)
        TEMPLATE(T, sparse_vec_set)(&B->rows[i], &M2->rows[i-M1->r], M2->c_off, ctx);
}

/* Split block matrix B = [M1 M2] into submatrices M1 and M2 */
/* M1->r and M2->r must equal B->r */
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_split_horizontal) (TEMPLATE(T, sparse_mat_t) M1, TEMPLATE(T, sparse_mat_t) M2, const TEMPLATE(T, sparse_mat_t) B, slong c, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    for(i=0; i<B->r; ++i) TEMPLATE(T, sparse_vec_split)(&M1->rows[i], &M2->rows[i], &B->rows[i], c, ctx);
}

/* Split block matix B = [M1^t M1^t]^t into submatrices M1 and M2 */
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_split_vertical) (TEMPLATE(T, sparse_mat_t) M1, TEMPLATE(T, sparse_mat_t) M2, const TEMPLATE(T, sparse_mat_t) B, slong r, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    r = FLINT_MIN(r, B->r);
    for(i=0; i<r; ++i) TEMPLATE(T, sparse_vec_set)(&M1->rows[i], &B->rows[i], B->c_off, ctx);
    for(i=r; i<B->r; ++i) TEMPLATE(T, sparse_vec_set)(&M2->rows[i-r], &B->rows[i], B->c_off, ctx);
}

/* Matrix permutation */
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_permute_cols)(TEMPLATE(T, sparse_mat_t) M, slong *Q, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < M->r; ++i) {
        if(!M->rows[i].nnz) continue;
        TEMPLATE(T, sparse_vec_permute_inds)(&M->rows[i], Q, ctx);
        qsort(M->rows[i].entries, M->rows[i].nnz, sizeof(*M->rows[i].entries), TEMPLATE(T, sparse_entry_cmp));
    }
}

FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_permute_rows)(TEMPLATE(T, sparse_mat_t) M, slong *P, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    TEMPLATE(T, sparse_vec_struct) *prows;
    prows = flint_calloc(M->r, sizeof(*prows));
    for (i = 0; i < M->r; ++i) prows[P[i]] = M->rows[i];
    memcpy(M->rows, prows, M->r*sizeof(*M->rows));
    flint_free(prows);
}

/* Random matrix generation */
FLINT_DLL void TEMPLATE(T, sparse_mat_randtest) (TEMPLATE(T, sparse_mat_t) M, flint_rand_t state, slong min_nnz, slong max_nnz, const TEMPLATE(T, ctx_t) ctx);
/*
FLINT_DLL void TEMPLATE(T, sparse_mat_randfull) (TEMPLATE(T, sparse_mat_t) M, flint_rand_t state, TEMPLATE(T, ctx_t) ctx);
FLINT_DLL int TEMPLATE(T, sparse_mat_randpermdiag)(TEMPLATE(T, sparse_mat_t) M, flint_rand_t state,
                 const TEMPLATE(T, struct) *diag, slong n);
FLINT_DLL void TEMPLATE(T, sparse_mat_randrank) (TEMPLATE(T, sparse_mat_t), flint_rand_t state, slong rank, TEMPLATE(T, ctx_t) ctx);
FLINT_DLL void TEMPLATE(T, sparse_mat_randops) (TEMPLATE(T, sparse_mat_t) M, slong count, flint_rand_t state, TEMPLATE(T, ctx_t) ctx);
FLINT_DLL void TEMPLATE(T, sparse_mat_randtril) (TEMPLATE(T, sparse_mat_t) M, flint_rand_t state, int unit, TEMPLATE(T, ctx_t) ctx);
FLINT_DLL void TEMPLATE(T, sparse_mat_randtriu) (TEMPLATE(T, sparse_mat_t) M, flint_rand_t state, int unit, TEMPLATE(T, ctx_t) ctx);
 */

FLINT_DLL void TEMPLATE(T, sparse_mat_print_pretty) (const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx);

FQ_SPARSE_MAT_TEMPLATES_INLINE
int TEMPLATE(T, sparse_mat_equal) (const TEMPLATE(T, sparse_mat_t) M1, const TEMPLATE(T, sparse_mat_t) M2, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    if (M1->r != M2->r) return 0;
    for (i = 0; i < M1->r; ++i)
        if (!TEMPLATE(T, sparse_vec_equal)(&M1->rows[i], &M2->rows[i], M1->c_off-M2->c_off, ctx)) return 0;
    return 1;
}

FQ_SPARSE_MAT_TEMPLATES_INLINE
int TEMPLATE(T, sparse_mat_is_zero) (const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < M->r; ++i) 
        if (!TEMPLATE(T, sparse_vec_is_zero)(&M->rows[i], ctx)) return 0;
    return 1;
}

/* Must have M->r == N->c and M->c == N->r */
FLINT_DLL void TEMPLATE(T, sparse_mat_transpose) (TEMPLATE(T, sparse_mat_t) N, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx);

/* Arithmetic */
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_neg) (TEMPLATE(T, sparse_mat_t) N, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < N->r; ++i) TEMPLATE(T, sparse_vec_neg)(&N->rows[i], &M->rows[i], ctx);
}

FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, TEMPLATE(sparse_mat_scalar_mul, T)) (TEMPLATE(T, sparse_mat_t) N, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, t) c, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < N->r; ++i) TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T))(&N->rows[i], &M->rows[i], c, ctx);    
}


FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_add) (TEMPLATE(T, sparse_mat_t) O, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, sparse_mat_t) N, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < O->r; ++i) TEMPLATE(T, sparse_vec_add)(&O->rows[i], &M->rows[i], &N->rows[i], ctx);
}

FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_sub) (TEMPLATE(T, sparse_mat_t) O, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, sparse_mat_t) N, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < O->r; ++i) TEMPLATE(T, sparse_vec_sub)(&O->rows[i], &M->rows[i], &N->rows[i], ctx);
}

FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, TEMPLATE(sparse_mat_scalar_addmul, T)) (TEMPLATE(T, sparse_mat_t) O, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, sparse_mat_t) N, const TEMPLATE(T, t) c, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < O->r; ++i) TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T))(&O->rows[i], &M->rows[i], &N->rows[i], c, ctx);
}

FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, TEMPLATE(sparse_mat_scalar_submul, T)) (TEMPLATE(T, sparse_mat_t) O, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, sparse_mat_t) N, const TEMPLATE(T, t) c, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < O->r; ++i) TEMPLATE(T, TEMPLATE(sparse_vec_scalar_submul, T))(&O->rows[i], &M->rows[i], &N->rows[i], c, ctx);
}

/* Matrix-vector and matrix-matrix multipliciation */
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_mul_vec) (TEMPLATE(T, struct) *y, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *x, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i;
    for (i = 0; i < M->r; ++i) TEMPLATE(T, sparse_vec_dot_dense)(&y[i], &M->rows[i], x, ctx);
}
FQ_SPARSE_MAT_TEMPLATES_INLINE
void TEMPLATE(T, sparse_mat_mul_mat) (TEMPLATE(T, mat_t) Y, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, mat_t) X, const TEMPLATE(T, ctx_t) ctx) 
{
    slong i, j;
    TEMPLATE(T, mat_zero) (Y, ctx);
    for (i = 0; i < M->r; ++i)
    {
        for (j = 0; j < M->rows[i].nnz; ++j)
        {
            TEMPLATE(T, sparse_entry_struct) *e = &M->rows[i].entries[j];
            _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T))(Y->rows[i], X->rows[e->ind], X->c, e->val, ctx);
        }
    }
}

FLINT_DLL
slong TEMPLATE(T, sparse_mat_inv) (TEMPLATE(T, sparse_mat_t) Mi, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx);

/* Decomposition/reduction */
FLINT_DLL
slong TEMPLATE(T, sparse_mat_lu)(slong *P, slong *Q, TEMPLATE(T, sparse_mat_t) L, TEMPLATE(T, sparse_mat_t) U, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
slong TEMPLATE(T, sparse_mat_rref) (TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx);

/* Solve Ax=b */
FLINT_DLL
int TEMPLATE(T, sparse_mat_solve_lanczos) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *b, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
int TEMPLATE(T, sparse_mat_solve_wiedemann) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *b, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
int TEMPLATE(T, sparse_mat_solve_lu) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *b, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
int TEMPLATE(T, sparse_mat_solve_rref) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *b, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
int TEMPLATE(T, sparse_mat_solve_block_wiedemann) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *b, slong block_size, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
int TEMPLATE(T, sparse_mat_solve_block_lanczos) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *b, slong block_size, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx);

/* Find single nullvector */
FLINT_DLL
int TEMPLATE(T, sparse_mat_nullvector_wiedemann) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
int TEMPLATE(T, sparse_mat_nullvector_lanczos) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
int TEMPLATE(T, sparse_mat_nullvector_block_wiedemann) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, slong block_size, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx); 

FLINT_DLL
int TEMPLATE(T, sparse_mat_nullvector_block_lanczos) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, slong block_size, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx); 

/* Note: this should take in uninitialized matrix X */
FLINT_DLL
slong TEMPLATE(T, sparse_mat_nullspace_rref) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
slong TEMPLATE(T, sparse_mat_nullspace_lu) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
slong TEMPLATE(T, sparse_mat_nullspace_lanczos) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, sparse_mat_t) M, flint_rand_t state, slong max_iters, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
slong TEMPLATE(T, sparse_mat_nullspace_wiedemann) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, sparse_mat_t) M, flint_rand_t state, slong max_iters, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
slong TEMPLATE(T, sparse_mat_nullspace_block_lanczos) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, sparse_mat_t) M, slong block_size, flint_rand_t state, slong max_iters, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL
slong TEMPLATE(T, sparse_mat_nullspace_block_wiedemann) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, sparse_mat_t) M, slong block_size, flint_rand_t state, slong max_iters, const TEMPLATE(T, ctx_t) ctx);

/* Nullspace */
/* NMOD_SPARSE_MAT_INLINE
slong TEMPLATE(T, sparse_mat_nullspace) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, sparse_mat_t) M, TEMPLATE(T, ctx_t) ctx);
 */
#ifdef __cplusplus
}
#endif

#endif

