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

#ifndef POLY_MAT_H
#define POLY_MAT_H

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Types *********************************************************************/

typedef struct
{
    fmpz_poly_struct * entries;
    long r;
    long c;
    fmpz_poly_struct ** rows;
} fmpz_poly_mat_struct;

typedef fmpz_poly_mat_struct fmpz_poly_mat_t[1];

#define fmpz_poly_mat_entry(mat,i,j) ((mat)->rows[(i)] + (j))

/* Memory management *********************************************************/

void fmpz_poly_mat_init(fmpz_poly_mat_t mat, long rows, long cols);

void fmpz_poly_mat_init_set(fmpz_poly_mat_t mat, const fmpz_poly_mat_t src);

void fmpz_poly_mat_swap(fmpz_poly_mat_t mat1, fmpz_poly_mat_t mat2);

void fmpz_poly_mat_set(fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2);

void fmpz_poly_mat_clear(fmpz_poly_mat_t mat);

/* Basic properties **********************************************************/

static __inline__ long
fmpz_poly_mat_nrows(const fmpz_poly_mat_t mat)
{
    return mat->r;
}

static __inline__ long
fmpz_poly_mat_ncols(const fmpz_poly_mat_t mat)
{
    return mat->c;
}

/* Comparison ****************************************************************/

int fmpz_poly_mat_equal(const fmpz_poly_mat_t mat1,
                        const fmpz_poly_mat_t mat2);

int fmpz_poly_mat_is_zero(const fmpz_poly_mat_t mat);

int fmpz_poly_mat_is_one(const fmpz_poly_mat_t mat);

static __inline__ int
fmpz_poly_mat_is_empty(const fmpz_poly_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

static __inline__ int
fmpz_poly_mat_is_square(const fmpz_poly_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Standard matrices *********************************************************/

void fmpz_poly_mat_zero(fmpz_poly_mat_t mat);

void fmpz_poly_mat_one(fmpz_poly_mat_t mat);

/* Random matrices ***********************************************************/

void fmpz_poly_mat_randtest(fmpz_poly_mat_t mat, flint_rand_t state,
                                long len, mp_bitcnt_t bits);

void fmpz_poly_mat_randtest_unsigned(fmpz_poly_mat_t mat, flint_rand_t state,
                             long len, mp_bitcnt_t bits);

void fmpz_poly_mat_randtest_sparse(fmpz_poly_mat_t A, flint_rand_t state,
                        long len, mp_bitcnt_t bits, float density);

/* Input and output **********************************************************/

void fmpz_poly_mat_print(const fmpz_poly_mat_t mat, const char * x);

/* Norms *********************************************************************/

long fmpz_poly_mat_max_bits(const fmpz_poly_mat_t A);

long fmpz_poly_mat_max_length(const fmpz_poly_mat_t A);

/* Transpose *****************************************************************/

void fmpz_poly_mat_transpose(fmpz_poly_mat_t B, const fmpz_poly_mat_t A);

/* Truncation ****************************************************************/

void fmpz_poly_mat_truncate(fmpz_poly_mat_t A, long len);

/* Scalar arithmetic *********************************************************/

void fmpz_poly_mat_scalar_mul_fmpz_poly(fmpz_poly_mat_t B,
                    const fmpz_poly_mat_t A, const fmpz_poly_t c);

void fmpz_poly_mat_scalar_mul_fmpz(fmpz_poly_mat_t B,
                    const fmpz_poly_mat_t A, const fmpz_t c);

/* Matrix arithmetic *********************************************************/

void fmpz_poly_mat_add(fmpz_poly_mat_t C, const fmpz_poly_mat_t A,
                                            const fmpz_poly_mat_t B);

void fmpz_poly_mat_sub(fmpz_poly_mat_t C, const fmpz_poly_mat_t A,
                                            const fmpz_poly_mat_t B);

void fmpz_poly_mat_neg(fmpz_poly_mat_t B, const fmpz_poly_mat_t A);

void fmpz_poly_mat_mul(fmpz_poly_mat_t C, const fmpz_poly_mat_t A,
                                            const fmpz_poly_mat_t B);

void fmpz_poly_mat_mul_classical(fmpz_poly_mat_t C, const fmpz_poly_mat_t A,
                                                const fmpz_poly_mat_t B);

void fmpz_poly_mat_mul_KS(fmpz_poly_mat_t C, const fmpz_poly_mat_t A,
                                            const fmpz_poly_mat_t B);

void fmpz_poly_mat_mullow(fmpz_poly_mat_t C, const fmpz_poly_mat_t A,
    const fmpz_poly_mat_t B, long len);

void fmpz_poly_mat_sqr(fmpz_poly_mat_t B, const fmpz_poly_mat_t A);

void fmpz_poly_mat_sqr_classical(fmpz_poly_mat_t B, const fmpz_poly_mat_t A);

void fmpz_poly_mat_sqr_KS(fmpz_poly_mat_t B, const fmpz_poly_mat_t A);

void fmpz_poly_mat_sqrlow(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, long len);

void fmpz_poly_mat_pow(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, ulong exp);

void
fmpz_poly_mat_pow_trunc(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, ulong exp,
                            long len);

void fmpz_poly_mat_prod(fmpz_poly_mat_t res,
                        fmpz_poly_mat_t * const factors, long n);

/* Evaluation ****************************************************************/

void fmpz_poly_mat_evaluate_fmpz(fmpz_mat_t B,
                        const fmpz_poly_mat_t A, const fmpz_t x);

/* Row reduction *************************************************************/

long fmpz_poly_mat_find_pivot_any(const fmpz_poly_mat_t mat,
                                    long start_row, long end_row, long c);

long fmpz_poly_mat_find_pivot_partial(const fmpz_poly_mat_t mat,
                                    long start_row, long end_row, long c);

long fmpz_poly_mat_fflu(fmpz_poly_mat_t B, fmpz_poly_t den, long * perm,
                            const fmpz_poly_mat_t A, int rank_check);

long fmpz_poly_mat_rref(fmpz_poly_mat_t B, fmpz_poly_t den,
                            const fmpz_poly_mat_t A);

/* Trace *********************************************************************/

void fmpz_poly_mat_trace(fmpz_poly_t trace, const fmpz_poly_mat_t mat);

/* Determinant and rank ******************************************************/

void fmpz_poly_mat_det(fmpz_poly_t det, const fmpz_poly_mat_t A);

void fmpz_poly_mat_det_fflu(fmpz_poly_t det, const fmpz_poly_mat_t A);

void fmpz_poly_mat_det_interpolate(fmpz_poly_t det, const fmpz_poly_mat_t A);

long fmpz_poly_mat_rank(const fmpz_poly_mat_t A);

/* Inverse *******************************************************************/

int fmpz_poly_mat_inv(fmpz_poly_mat_t Ainv, fmpz_poly_t den,
    const fmpz_poly_mat_t A);

/* Nullspace *****************************************************************/

long fmpz_poly_mat_nullspace(fmpz_poly_mat_t res, const fmpz_poly_mat_t mat);

/* Solving *******************************************************************/

int fmpz_poly_mat_solve(fmpz_poly_mat_t X, fmpz_poly_t den,
                    const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

int fmpz_poly_mat_solve_fflu(fmpz_poly_mat_t X, fmpz_poly_t den,
                            const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

void fmpz_poly_mat_solve_fflu_precomp(fmpz_poly_mat_t X,
                    const long * perm,
                    const fmpz_poly_mat_t FFLU, const fmpz_poly_mat_t B);

#ifdef __cplusplus
}
#endif

#endif
