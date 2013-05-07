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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#ifndef NMOD_POLY_MAT_H
#define NMOD_POLY_MAT_H

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Types *********************************************************************/

typedef struct
{
    nmod_poly_struct * entries;
    long r;
    long c;
    nmod_poly_struct ** rows;
    mp_limb_t modulus;
} nmod_poly_mat_struct;

typedef nmod_poly_mat_struct nmod_poly_mat_t[1];

#define nmod_poly_mat_entry(mat,i,j) ((mat)->rows[(i)] + (j))

/* Memory management *********************************************************/

void nmod_poly_mat_init(nmod_poly_mat_t mat, long rows, long cols, mp_limb_t n);

void nmod_poly_mat_init_set(nmod_poly_mat_t mat, const nmod_poly_mat_t src);

void nmod_poly_mat_swap(nmod_poly_mat_t mat1, nmod_poly_mat_t mat2);

void nmod_poly_mat_set(nmod_poly_mat_t mat1, const nmod_poly_mat_t mat2);

void nmod_poly_mat_clear(nmod_poly_mat_t mat);

/* Basic properties **********************************************************/

static __inline__ long
nmod_poly_mat_nrows(const nmod_poly_mat_t mat)
{
    return mat->r;
}

static __inline__ long
nmod_poly_mat_ncols(const nmod_poly_mat_t mat)
{
    return mat->c;
}

static __inline__ mp_limb_t
nmod_poly_mat_modulus(const nmod_poly_mat_t mat)
{
    return mat->modulus;
}

/* Comparison ****************************************************************/

int nmod_poly_mat_equal(const nmod_poly_mat_t mat1,
                        const nmod_poly_mat_t mat2);

int nmod_poly_mat_is_zero(const nmod_poly_mat_t mat);

int nmod_poly_mat_is_one(const nmod_poly_mat_t mat);

static __inline__ int
nmod_poly_mat_is_empty(const nmod_poly_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

static __inline__ int
nmod_poly_mat_is_square(const nmod_poly_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Standard matrices *********************************************************/

void nmod_poly_mat_zero(nmod_poly_mat_t mat);

void nmod_poly_mat_one(nmod_poly_mat_t mat);

/* Random matrices ***********************************************************/

void nmod_poly_mat_randtest(nmod_poly_mat_t mat, flint_rand_t state,
                                long len);

void nmod_poly_mat_randtest_sparse(nmod_poly_mat_t A, flint_rand_t state,
                        long len, float density);

/* Input and output **********************************************************/

void nmod_poly_mat_print(const nmod_poly_mat_t mat, const char * x);

/* Norms *********************************************************************/

long nmod_poly_mat_max_length(const nmod_poly_mat_t A);

/* Scalar arithmetic *********************************************************/

void nmod_poly_mat_scalar_mul_nmod_poly(nmod_poly_mat_t B,
                    const nmod_poly_mat_t A, const nmod_poly_t c);

void nmod_poly_mat_scalar_mul_nmod(nmod_poly_mat_t B,
                    const nmod_poly_mat_t A, mp_limb_t c);

/* Matrix arithmetic *********************************************************/

void nmod_poly_mat_add(nmod_poly_mat_t C, const nmod_poly_mat_t A,
                                            const nmod_poly_mat_t B);

void nmod_poly_mat_sub(nmod_poly_mat_t C, const nmod_poly_mat_t A,
                                            const nmod_poly_mat_t B);

void nmod_poly_mat_neg(nmod_poly_mat_t B, const nmod_poly_mat_t A);

void nmod_poly_mat_mul(nmod_poly_mat_t C, const nmod_poly_mat_t A,
                                            const nmod_poly_mat_t B);

void nmod_poly_mat_mul_interpolate(nmod_poly_mat_t C, const nmod_poly_mat_t A,
    const nmod_poly_mat_t B);

void nmod_poly_mat_mul_classical(nmod_poly_mat_t C, const nmod_poly_mat_t A,
                                            const nmod_poly_mat_t B);

void nmod_poly_mat_mul_KS(nmod_poly_mat_t C, const nmod_poly_mat_t A,
    const nmod_poly_mat_t B);

void nmod_poly_mat_sqr(nmod_poly_mat_t B, const nmod_poly_mat_t A);

void nmod_poly_mat_sqr_classical(nmod_poly_mat_t B, const nmod_poly_mat_t A);

void nmod_poly_mat_sqr_KS(nmod_poly_mat_t B, const nmod_poly_mat_t A);

void nmod_poly_mat_sqr_interpolate(nmod_poly_mat_t B, const nmod_poly_mat_t A);

void nmod_poly_mat_pow(nmod_poly_mat_t B, const nmod_poly_mat_t A, ulong exp);

/* Evaluation ****************************************************************/

void nmod_poly_mat_evaluate_nmod(nmod_mat_t B, const nmod_poly_mat_t A, mp_limb_t x);

/* Row reduction *************************************************************/

long nmod_poly_mat_find_pivot_any(const nmod_poly_mat_t mat,
                                    long start_row, long end_row, long c);

long nmod_poly_mat_find_pivot_partial(const nmod_poly_mat_t mat,
                                    long start_row, long end_row, long c);

long nmod_poly_mat_fflu(nmod_poly_mat_t B, nmod_poly_t den, long * perm,
                            const nmod_poly_mat_t A, int rank_check);

long nmod_poly_mat_rref(nmod_poly_mat_t B, nmod_poly_t den,
                            const nmod_poly_mat_t A);

/* Trace *********************************************************************/

void nmod_poly_mat_trace(nmod_poly_t trace, const nmod_poly_mat_t mat);

/* Determinant and rank ******************************************************/

void nmod_poly_mat_det(nmod_poly_t det, const nmod_poly_mat_t A);

void nmod_poly_mat_det_fflu(nmod_poly_t det, const nmod_poly_mat_t A);

void nmod_poly_mat_det_interpolate(nmod_poly_t det, const nmod_poly_mat_t A);

long nmod_poly_mat_rank(const nmod_poly_mat_t A);

/* Inverse *******************************************************************/

int nmod_poly_mat_inv(nmod_poly_mat_t Ainv, nmod_poly_t den,
    const nmod_poly_mat_t A);

/* Nullspace *****************************************************************/

long nmod_poly_mat_nullspace(nmod_poly_mat_t res, const nmod_poly_mat_t mat);

/* Solving *******************************************************************/

int nmod_poly_mat_solve(nmod_poly_mat_t X, nmod_poly_t den,
                    const nmod_poly_mat_t A, const nmod_poly_mat_t B);

int nmod_poly_mat_solve_fflu(nmod_poly_mat_t X, nmod_poly_t den,
                            const nmod_poly_mat_t A, const nmod_poly_mat_t B);

void nmod_poly_mat_solve_fflu_precomp(nmod_poly_mat_t X,
                    const long * perm,
                    const nmod_poly_mat_t FFLU, const nmod_poly_mat_t B);

#ifdef __cplusplus
}
#endif

#endif
