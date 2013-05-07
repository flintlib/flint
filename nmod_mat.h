/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 William Hart
    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/

#ifndef NMOD_MAT_H
#define NMOD_MAT_H

#undef ulong /* interferes with system includes */
#include <stdlib.h>
#define ulong unsigned long

#include <gmp.h>
#include "flint.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    mp_limb_t * entries;
    long r;
    long c;
    mp_limb_t ** rows;
    nmod_t mod;
}
nmod_mat_struct;

/* nmod_mat_t allows reference-like semantics for nmod_mat_struct */
typedef nmod_mat_struct nmod_mat_t[1];

#define nmod_mat_entry(mat,i,j) ((mat)->rows[(i)][(j)])
#define nmod_mat_nrows(mat) ((mat)->r)
#define nmod_mat_ncols(mat) ((mat)->c)

static __inline__
void
_nmod_mat_set_mod(nmod_mat_t mat, mp_limb_t n)
{
    mat->mod.n = n;
    mat->mod.ninv = n_preinvert_limb(n);
    count_leading_zeros(mat->mod.norm, n);
}

/* Memory management */
void nmod_mat_init(nmod_mat_t mat, long rows, long cols, mp_limb_t n);
void nmod_mat_init_set(nmod_mat_t mat, const nmod_mat_t src);
void nmod_mat_clear(nmod_mat_t mat);

void nmod_mat_window_init(nmod_mat_t window, const nmod_mat_t mat, long r1, long c1, long r2, long c2);
void nmod_mat_window_clear(nmod_mat_t window);

/* Random matrix generation */
void nmod_mat_randtest(nmod_mat_t mat, flint_rand_t state);
void nmod_mat_randfull(nmod_mat_t mat, flint_rand_t state);
int nmod_mat_randpermdiag(nmod_mat_t mat, flint_rand_t state, 
                 mp_srcptr diag, long n);
void nmod_mat_randrank(nmod_mat_t, flint_rand_t state, long rank);
void nmod_mat_randops(nmod_mat_t mat, long count, flint_rand_t state);
void nmod_mat_randtril(nmod_mat_t mat, flint_rand_t state, int unit);
void nmod_mat_randtriu(nmod_mat_t mat, flint_rand_t state, int unit);


void nmod_mat_print_pretty(const nmod_mat_t mat);

int nmod_mat_equal(const nmod_mat_t mat1, const nmod_mat_t mat2);

void nmod_mat_zero(nmod_mat_t mat);

int nmod_mat_is_zero(const nmod_mat_t mat);

static __inline__ int
nmod_mat_is_empty(const nmod_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

static __inline__ int
nmod_mat_is_square(const nmod_mat_t mat)
{
    return (mat->r == mat->c);
}


void nmod_mat_set(nmod_mat_t B, const nmod_mat_t A);
void nmod_mat_transpose(nmod_mat_t B, const nmod_mat_t A);

/* Addition and subtraction */

void nmod_mat_add(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void nmod_mat_sub(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void nmod_mat_neg(nmod_mat_t B, const nmod_mat_t A);

/* Matrix-scalar arithmetic */

void nmod_mat_scalar_mul(nmod_mat_t B, const nmod_mat_t A, mp_limb_t c);

/* Matrix multiplication */

void nmod_mat_mul(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void nmod_mat_mul_classical(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void nmod_mat_mul_strassen(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);

void
_nmod_mat_mul_classical(nmod_mat_t D, const nmod_mat_t C,
                                const nmod_mat_t A, const nmod_mat_t B, int op);

void nmod_mat_addmul(nmod_mat_t D, const nmod_mat_t C,
                                const nmod_mat_t A, const nmod_mat_t B);

void nmod_mat_submul(nmod_mat_t D, const nmod_mat_t C,
                                const nmod_mat_t A, const nmod_mat_t B);

/* Trace */

mp_limb_t nmod_mat_trace(const nmod_mat_t mat);

/* Determinant */

mp_limb_t _nmod_mat_det(nmod_mat_t A);
mp_limb_t nmod_mat_det(const nmod_mat_t A);

/* Rank */

long nmod_mat_rank(const nmod_mat_t A);

/* Inverse */

int nmod_mat_inv(nmod_mat_t B, const nmod_mat_t A);

/* Triangular solving */

void nmod_mat_solve_tril(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);
void nmod_mat_solve_tril_recursive(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);
void nmod_mat_solve_tril_classical(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);

void nmod_mat_solve_triu(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit);
void nmod_mat_solve_triu_recursive(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit);
void nmod_mat_solve_triu_classical(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit);

/* LU decomposition */

long nmod_mat_lu(long * P, nmod_mat_t A, int rank_check);
long nmod_mat_lu_classical(long * P, nmod_mat_t A, int rank_check);
long nmod_mat_lu_recursive(long * P, nmod_mat_t A, int rank_check);

/* Nonsingular solving */

int nmod_mat_solve(nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B);
int nmod_mat_solve_vec(mp_ptr x, const nmod_mat_t A, mp_srcptr b);

/* Reduced row echelon form */

long nmod_mat_rref(nmod_mat_t A);

/* Nullspace */

long nmod_mat_nullspace(nmod_mat_t X, const nmod_mat_t A);


/* Tuning parameters *********************************************************/

/* Size at which pre-transposing becomes faster in classical multiplication */
#define NMOD_MAT_MUL_TRANSPOSE_CUTOFF 20

/* Strassen multiplication */
#define NMOD_MAT_MUL_STRASSEN_CUTOFF 256

/* Cutoff between classical and recursive triangular solving */
#define NMOD_MAT_SOLVE_TRI_ROWS_CUTOFF 64
#define NMOD_MAT_SOLVE_TRI_COLS_CUTOFF 64

/* Cutoff between classical and recursive LU decomposition */
#define NMOD_MAT_LU_RECURSIVE_CUTOFF 4

/*
   Suggested initial modulus size for multimodular algorithms. This should
   be chosen so that we get the most number of bits per cycle
   in matrix multiplication. On x86-64 it appears to be optimal to use
   moduli giving nlimbs = 2. This should hold both in the classical
   range and in Strassen blocks.
 */
#define NMOD_MAT_OPTIMAL_MODULUS_BITS (FLINT_BITS-5)

#ifdef __cplusplus
}
#endif

#endif

