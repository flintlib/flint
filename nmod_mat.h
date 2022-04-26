/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010,2011 Fredrik Johansson
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MAT_H
#define NMOD_MAT_H

#ifdef NMOD_MAT_INLINES_C
#define NMOD_MAT_INLINE FLINT_DLL
#else
#define NMOD_MAT_INLINE static __inline__
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
#include "fmpz.h"
#include "thread_support.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    mp_limb_t * entries;
    slong r;
    slong c;
    mp_limb_t ** rows;
    nmod_t mod;
}
nmod_mat_struct;

/* nmod_mat_t allows reference-like semantics for nmod_mat_struct */
typedef nmod_mat_struct nmod_mat_t[1];

#define nmod_mat_entry(mat,i,j) ((mat)->rows[(i)][(j)])

NMOD_MAT_INLINE
mp_limb_t nmod_mat_get_entry(const nmod_mat_t mat, slong i, slong j)
{
   return mat->rows[i][j];
}

NMOD_MAT_INLINE
mp_limb_t * nmod_mat_entry_ptr(const nmod_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

/* See inlines.c */

NMOD_MAT_INLINE
slong nmod_mat_nrows(const nmod_mat_t mat)
{
   return mat->r;
}

NMOD_MAT_INLINE
slong nmod_mat_ncols(const nmod_mat_t mat)
{
   return mat->c;
}

NMOD_MAT_INLINE
void _nmod_mat_set_mod(nmod_mat_t mat, mp_limb_t n)
{
    mat->mod.n = n;
    count_leading_zeros(mat->mod.norm, n);
    invert_limb(mat->mod.ninv, n << mat->mod.norm);
}

/* Memory management */
FLINT_DLL void nmod_mat_init(nmod_mat_t mat, slong rows, slong cols, mp_limb_t n);
FLINT_DLL void nmod_mat_init_set(nmod_mat_t mat, const nmod_mat_t src);
FLINT_DLL void nmod_mat_clear(nmod_mat_t mat);
FLINT_DLL void nmod_mat_one(nmod_mat_t mat);

FLINT_DLL void nmod_mat_swap(nmod_mat_t mat1, nmod_mat_t mat2);

NMOD_MAT_INLINE void
nmod_mat_swap_entrywise(nmod_mat_t mat1, nmod_mat_t mat2)
{
    slong i, j;
    for (i = 0; i < nmod_mat_nrows(mat1); i++)
    {
       mp_limb_t * row1 = mat1->rows[i];
       mp_limb_t * row2 = mat2->rows[i];
       for (j = 0; j < nmod_mat_ncols(mat1); j++)
          MP_LIMB_SWAP(row1[j], row2[j]);
    }
}

/* Windows and concatenation */

FLINT_DLL void nmod_mat_window_init(nmod_mat_t window, const nmod_mat_t mat, slong r1, slong c1, slong r2, slong c2);
FLINT_DLL void nmod_mat_window_clear(nmod_mat_t window);

FLINT_DLL void nmod_mat_concat_horizontal(nmod_mat_t res,
                           const nmod_mat_t mat1,  const nmod_mat_t mat2);
FLINT_DLL void nmod_mat_concat_vertical(nmod_mat_t res,
                           const nmod_mat_t mat1,  const nmod_mat_t mat2);

/* Random matrix generation */
FLINT_DLL void nmod_mat_randtest(nmod_mat_t mat, flint_rand_t state);
FLINT_DLL void nmod_mat_randfull(nmod_mat_t mat, flint_rand_t state);
FLINT_DLL int nmod_mat_randpermdiag(nmod_mat_t mat, flint_rand_t state,
                 mp_srcptr diag, slong n);
FLINT_DLL void nmod_mat_randrank(nmod_mat_t, flint_rand_t state, slong rank);
FLINT_DLL void nmod_mat_randops(nmod_mat_t mat, slong count, flint_rand_t state);
FLINT_DLL void nmod_mat_randtril(nmod_mat_t mat, flint_rand_t state, int unit);
FLINT_DLL void nmod_mat_randtriu(nmod_mat_t mat, flint_rand_t state, int unit);


FLINT_DLL int nmod_mat_fprint_pretty(FILE* file, const nmod_mat_t mat);

NMOD_MAT_INLINE void nmod_mat_print_pretty(const nmod_mat_t mat)
{
    nmod_mat_fprint_pretty(stdout, mat);
}

NMOD_MAT_INLINE int nmod_mat_print(const nmod_mat_t mat)
{
    return nmod_mat_fprint_pretty(stdout, mat);
}

NMOD_MAT_INLINE int nmod_mat_fprint(FILE* f, const nmod_mat_t mat)
{
    return nmod_mat_fprint_pretty(f, mat);
}

FLINT_DLL int nmod_mat_equal(const nmod_mat_t mat1, const nmod_mat_t mat2);

FLINT_DLL void nmod_mat_zero(nmod_mat_t mat);

FLINT_DLL int nmod_mat_is_zero(const nmod_mat_t mat);

FLINT_DLL int nmod_mat_is_one(const nmod_mat_t mat);

NMOD_MAT_INLINE
int nmod_mat_is_zero_row(const nmod_mat_t mat, slong i)
{
    return _nmod_vec_is_zero(mat->rows[i], mat->c);
}

NMOD_MAT_INLINE
int nmod_mat_is_empty(const nmod_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

NMOD_MAT_INLINE
int nmod_mat_is_square(const nmod_mat_t mat)
{
    return (mat->r == mat->c);
}


FLINT_DLL void nmod_mat_set(nmod_mat_t B, const nmod_mat_t A);
FLINT_DLL void nmod_mat_transpose(nmod_mat_t B, const nmod_mat_t A);

/* Addition and subtraction */

FLINT_DLL void nmod_mat_add(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
FLINT_DLL void nmod_mat_sub(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
FLINT_DLL void nmod_mat_neg(nmod_mat_t B, const nmod_mat_t A);

/* Matrix-scalar arithmetic */

FLINT_DLL void nmod_mat_scalar_mul(nmod_mat_t B, const nmod_mat_t A, mp_limb_t c);
FLINT_DLL void nmod_mat_scalar_addmul_ui(nmod_mat_t dest,
                       const nmod_mat_t X, const nmod_mat_t Y, const mp_limb_t b);


NMOD_MAT_INLINE
void nmod_mat_scalar_mul_fmpz(nmod_mat_t res, const nmod_mat_t M, const fmpz_t c)
{
    fmpz_t d;
    fmpz_init(d);
    fmpz_mod_ui(d, c, res->mod.n);
    nmod_mat_scalar_mul(res, M, fmpz_get_ui(d));
    fmpz_clear(d);
}

/* Matrix multiplication */

FLINT_DLL void nmod_mat_mul(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);

FLINT_DLL int nmod_mat_mul_blas(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);

FLINT_DLL void nmod_mat_mul_classical(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);

FLINT_DLL void
_nmod_mat_mul_classical_threaded_pool_op(nmod_mat_t D, const nmod_mat_t C,
		            const nmod_mat_t A, const nmod_mat_t B, int op,
			      thread_pool_handle * threads, slong num_threads);

FLINT_DLL void nmod_mat_mul_classical_threaded(nmod_mat_t C,
		                       const nmod_mat_t A, const nmod_mat_t B);
FLINT_DLL void nmod_mat_mul_strassen(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);

FLINT_DLL void _nmod_mat_mul_classical_op(nmod_mat_t D, const nmod_mat_t C,
                                const nmod_mat_t A, const nmod_mat_t B, int op);

FLINT_DLL void nmod_mat_addmul(nmod_mat_t D, const nmod_mat_t C,
                                const nmod_mat_t A, const nmod_mat_t B);

FLINT_DLL void nmod_mat_submul(nmod_mat_t D, const nmod_mat_t C,
                                const nmod_mat_t A, const nmod_mat_t B);

FLINT_DLL void nmod_mat_mul_nmod_vec(mp_limb_t * c, const nmod_mat_t A,
                                              const mp_limb_t * b, slong blen);

FLINT_DLL void nmod_mat_mul_nmod_vec_ptr(mp_limb_t * const * c,
                  const nmod_mat_t A, const mp_limb_t * const * b, slong blen);

FLINT_DLL void nmod_mat_nmod_vec_mul(mp_limb_t * c, const mp_limb_t * a,
                                               slong alen, const nmod_mat_t B);

FLINT_DLL void nmod_mat_nmod_vec_mul_ptr(mp_limb_t * const * c,
                  const mp_limb_t * const * a, slong alen, const nmod_mat_t B);

/* Exponent */

FLINT_DLL void _nmod_mat_pow(nmod_mat_t dest, const nmod_mat_t mat, ulong pow);
FLINT_DLL void nmod_mat_pow(nmod_mat_t dest, const nmod_mat_t mat, ulong pow);

/* Trace */

FLINT_DLL mp_limb_t nmod_mat_trace(const nmod_mat_t mat);

/* Determinant */

FLINT_DLL mp_limb_t _nmod_mat_det(nmod_mat_t A);
FLINT_DLL mp_limb_t nmod_mat_det(const nmod_mat_t A);

FLINT_DLL mp_limb_t _nmod_mat_det_howell(nmod_mat_t A);
FLINT_DLL mp_limb_t nmod_mat_det_howell(const nmod_mat_t A);

/* Rank */

FLINT_DLL slong nmod_mat_rank(const nmod_mat_t A);

/* Inverse */

FLINT_DLL int nmod_mat_inv(nmod_mat_t B, const nmod_mat_t A);

/* Permutations */

NMOD_MAT_INLINE
void nmod_mat_swap_rows(nmod_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s && !nmod_mat_is_empty(mat))
    {
        mp_limb_t * u;
        slong t;

        if (perm)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u;
    }
}

NMOD_MAT_INLINE
void nmod_mat_invert_rows(nmod_mat_t mat, slong * perm)
{
    slong i;

    for (i = 0; i < mat->r/2; i++)
        nmod_mat_swap_rows(mat, perm, i, mat->r - i - 1);
}

NMOD_MAT_INLINE
void nmod_mat_swap_cols(nmod_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s && !nmod_mat_is_empty(mat))
    {
        slong t;

        if (perm)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        for (t = 0; t < mat->r; t++)
        {
            mp_limb_t c = mat->rows[t][r];
            mat->rows[t][r] = mat->rows[t][s];
            mat->rows[t][s] = c;
        }
    }
}

NMOD_MAT_INLINE
void nmod_mat_invert_cols(nmod_mat_t mat, slong * perm)
{
    if (!(nmod_mat_is_empty(mat)))
    {
        slong t;
        slong i;
        slong c = mat->c;
        slong k = mat->c/2;
        mp_limb_t e;

        if (perm)
        {
            for (i =0; i < k; i++)
            {
                t = perm[i];
                perm[i] = perm[c - i];
                perm[c - i] = t;
            }
        }

        for (t = 0; t < mat->r; t++)
        {
            for (i = 0; i < k; i++)
            {
                e = mat->rows[t][i];
                mat->rows[t][i] = mat->rows[t][c - i - 1];
                mat->rows[t][c - i - 1] = e;
            }
        }
    }
}

FLINT_DLL void nmod_mat_apply_permutation(nmod_mat_t A, slong * P, slong n);

/* Triangular solving */

FLINT_DLL void nmod_mat_solve_tril(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);
FLINT_DLL void nmod_mat_solve_tril_recursive(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);
FLINT_DLL void nmod_mat_solve_tril_classical(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);

FLINT_DLL void nmod_mat_solve_triu(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit);
FLINT_DLL void nmod_mat_solve_triu_recursive(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit);
FLINT_DLL void nmod_mat_solve_triu_classical(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit);

/* LU decomposition */

FLINT_DLL slong nmod_mat_lu(slong * P, nmod_mat_t A, int rank_check);
FLINT_DLL slong nmod_mat_lu_classical(slong * P, nmod_mat_t A, int rank_check);
FLINT_DLL slong nmod_mat_lu_classical_delayed(slong * P, nmod_mat_t A, int rank_check);
FLINT_DLL slong nmod_mat_lu_recursive(slong * P, nmod_mat_t A, int rank_check);

/* Nonsingular solving */

FLINT_DLL int nmod_mat_solve(nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B);
FLINT_DLL int nmod_mat_solve_vec(mp_ptr x, const nmod_mat_t A, mp_srcptr b);

/* Solving */

FLINT_DLL int nmod_mat_can_solve_inner(slong * rank, slong * prm, slong * piv,
                          nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B);
    
FLINT_DLL int nmod_mat_can_solve(nmod_mat_t X, const nmod_mat_t A,
                                                            const nmod_mat_t B);

/* Reduced row echelon form */

FLINT_DLL slong nmod_mat_rref(nmod_mat_t A);
FLINT_DLL slong _nmod_mat_rref(nmod_mat_t A, slong * pivots_nonpivots, slong * P);
FLINT_DLL slong nmod_mat_rref_classical(nmod_mat_t A);
FLINT_DLL slong _nmod_mat_rref_classical(nmod_mat_t A, slong * pivots_nonpivots);
FLINT_DLL slong nmod_mat_rref_storjohann(nmod_mat_t A);
FLINT_DLL slong _nmod_mat_rref_storjohann(nmod_mat_t A, slong * pivots_nonpivots);

FLINT_DLL slong nmod_mat_reduce_row(nmod_mat_t M, slong * P, slong * L, slong m);

/* Nullspace */

FLINT_DLL slong nmod_mat_nullspace(nmod_mat_t X, const nmod_mat_t A);

/* Howell form */

FLINT_DLL void nmod_mat_strong_echelon_form(nmod_mat_t A);

FLINT_DLL slong nmod_mat_howell_form(nmod_mat_t A);

/* Transforms */

FLINT_DLL void nmod_mat_similarity(nmod_mat_t M, slong r, ulong d);

/* Characteristic polynomial and minimal polynomial */

/* The following prototype actually lives in nmod_poly.h
 * 
 * FLINT_DLL void nmod_mat_charpoly_danilevsky(nmod_poly_t p, const nmod_mat_t M);
 *
 * FLINT_DLL void nmod_mat_minpoly(nmod_poly_t p, const nmod_mat_t M);
*/

/* Tuning parameters *********************************************************/

/* Size at which pre-transposing becomes faster in classical multiplication */
#define NMOD_MAT_MUL_TRANSPOSE_CUTOFF 20

/* Cutoff between classical and recursive triangular solving */
#define NMOD_MAT_SOLVE_TRI_ROWS_CUTOFF 64
#define NMOD_MAT_SOLVE_TRI_COLS_CUTOFF 64

/*
   Suggested initial modulus size for multimodular algorithms. This should
   be chosen so that we get the most number of bits per cycle
   in matrix multiplication. On x86-64 it appears to be optimal to use
   moduli giving nlimbs = 2. This should hold both in the classical
   range and in Strassen blocks.
 */
#define NMOD_MAT_OPTIMAL_MODULUS_BITS (FLINT_BITS-5)

/* Inlines *******************************************************************/

FLINT_DLL void nmod_mat_set_entry(nmod_mat_t mat, slong i, slong j, mp_limb_t x);

#ifdef __cplusplus
}
#endif

#endif

