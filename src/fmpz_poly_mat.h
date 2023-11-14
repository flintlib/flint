/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef POLY_MAT_H
#define POLY_MAT_H

#ifdef FMPZ_POLY_MAT_INLINES_C
#define FMPZ_POLY_MAT_INLINE
#else
#define FMPZ_POLY_MAT_INLINE static inline
#endif

#include "fmpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Types *********************************************************************/

FMPZ_POLY_MAT_INLINE
fmpz_poly_struct * fmpz_poly_mat_entry(const fmpz_poly_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

FMPZ_POLY_MAT_INLINE
slong fmpz_poly_mat_nrows(const fmpz_poly_mat_t mat)
{
    return mat->r;
}

FMPZ_POLY_MAT_INLINE
slong fmpz_poly_mat_ncols(const fmpz_poly_mat_t mat)
{
    return mat->c;
}

/* Memory management *********************************************************/

void fmpz_poly_mat_init(fmpz_poly_mat_t mat, slong rows, slong cols);

void fmpz_poly_mat_init_set(fmpz_poly_mat_t mat, const fmpz_poly_mat_t src);

void fmpz_poly_mat_swap(fmpz_poly_mat_t mat1, fmpz_poly_mat_t mat2);

FMPZ_POLY_MAT_INLINE void
fmpz_poly_mat_swap_entrywise(fmpz_poly_mat_t mat1, fmpz_poly_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < fmpz_poly_mat_nrows(mat1); i++)
        for (j = 0; j < fmpz_poly_mat_ncols(mat1); j++)
            fmpz_poly_swap(fmpz_poly_mat_entry(mat2, i, j), fmpz_poly_mat_entry(mat1, i, j));
}

void fmpz_poly_mat_set(fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2);

void fmpz_poly_mat_clear(fmpz_poly_mat_t mat);

/* Comparison ****************************************************************/

int fmpz_poly_mat_equal(const fmpz_poly_mat_t mat1,
                        const fmpz_poly_mat_t mat2);

int fmpz_poly_mat_is_zero(const fmpz_poly_mat_t mat);

int fmpz_poly_mat_is_one(const fmpz_poly_mat_t mat);

FMPZ_POLY_MAT_INLINE
int fmpz_poly_mat_is_empty(const fmpz_poly_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

FMPZ_POLY_MAT_INLINE
int fmpz_poly_mat_is_square(const fmpz_poly_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Standard matrices *********************************************************/

void fmpz_poly_mat_zero(fmpz_poly_mat_t mat);

void fmpz_poly_mat_one(fmpz_poly_mat_t mat);

/* Random matrices ***********************************************************/

void fmpz_poly_mat_randtest(fmpz_poly_mat_t mat, flint_rand_t state,
                                slong len, flint_bitcnt_t bits);

void fmpz_poly_mat_randtest_unsigned(fmpz_poly_mat_t mat, flint_rand_t state,
                             slong len, flint_bitcnt_t bits);

void fmpz_poly_mat_randtest_sparse(fmpz_poly_mat_t A, flint_rand_t state,
                        slong len, flint_bitcnt_t bits, float density);

/* Windows and concatenation */

void fmpz_poly_mat_window_init(fmpz_poly_mat_t window, const fmpz_poly_mat_t mat, slong r1,
    slong c1, slong r2, slong c2);

void fmpz_poly_mat_window_clear(fmpz_poly_mat_t window);

void fmpz_poly_mat_concat_horizontal(fmpz_poly_mat_t res,
                           const fmpz_poly_mat_t mat1,  const fmpz_poly_mat_t mat2);

void fmpz_poly_mat_concat_vertical(fmpz_poly_mat_t res,
                           const fmpz_poly_mat_t mat1,  const fmpz_poly_mat_t mat2);

/* Input and output **********************************************************/

void fmpz_poly_mat_print(const fmpz_poly_mat_t mat, const char * x);

/* Norms *********************************************************************/

slong fmpz_poly_mat_max_bits(const fmpz_poly_mat_t A);

slong fmpz_poly_mat_max_length(const fmpz_poly_mat_t A);

/* Transpose *****************************************************************/

void fmpz_poly_mat_transpose(fmpz_poly_mat_t B, const fmpz_poly_mat_t A);

/* Truncation ****************************************************************/

void fmpz_poly_mat_truncate(fmpz_poly_mat_t A, slong len);

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
    const fmpz_poly_mat_t B, slong len);

void fmpz_poly_mat_sqr(fmpz_poly_mat_t B, const fmpz_poly_mat_t A);

void fmpz_poly_mat_sqr_classical(fmpz_poly_mat_t B, const fmpz_poly_mat_t A);

void fmpz_poly_mat_sqr_KS(fmpz_poly_mat_t B, const fmpz_poly_mat_t A);

void fmpz_poly_mat_sqrlow(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, slong len);

void fmpz_poly_mat_pow(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, ulong exp);

void fmpz_poly_mat_pow_trunc(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, ulong exp,
                            slong len);

void fmpz_poly_mat_prod(fmpz_poly_mat_t res,
                        fmpz_poly_mat_t * const factors, slong n);

/* Evaluation ****************************************************************/

void fmpz_poly_mat_evaluate_fmpz(fmpz_mat_t B,
                        const fmpz_poly_mat_t A, const fmpz_t x);

/* Row reduction *************************************************************/

slong fmpz_poly_mat_find_pivot_any(const fmpz_poly_mat_t mat,
                                    slong start_row, slong end_row, slong c);

slong fmpz_poly_mat_find_pivot_partial(const fmpz_poly_mat_t mat,
                                    slong start_row, slong end_row, slong c);

slong fmpz_poly_mat_fflu(fmpz_poly_mat_t B, fmpz_poly_t den, slong * perm,
                            const fmpz_poly_mat_t A, int rank_check);

slong fmpz_poly_mat_rref(fmpz_poly_mat_t B, fmpz_poly_t den,
                            const fmpz_poly_mat_t A);

/* Trace *********************************************************************/

void fmpz_poly_mat_trace(fmpz_poly_t trace, const fmpz_poly_mat_t mat);

/* Determinant and rank ******************************************************/

void fmpz_poly_mat_det(fmpz_poly_t det, const fmpz_poly_mat_t A);

void fmpz_poly_mat_det_fflu(fmpz_poly_t det, const fmpz_poly_mat_t A);

void fmpz_poly_mat_det_interpolate(fmpz_poly_t det, const fmpz_poly_mat_t A);

slong fmpz_poly_mat_rank(const fmpz_poly_mat_t A);

/* Inverse *******************************************************************/

int fmpz_poly_mat_inv(fmpz_poly_mat_t Ainv, fmpz_poly_t den,
    const fmpz_poly_mat_t A);

/* Nullspace *****************************************************************/

slong fmpz_poly_mat_nullspace(fmpz_poly_mat_t res, const fmpz_poly_mat_t mat);

/* Solving *******************************************************************/

int fmpz_poly_mat_solve(fmpz_poly_mat_t X, fmpz_poly_t den,
                    const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

int fmpz_poly_mat_solve_fflu(fmpz_poly_mat_t X, fmpz_poly_t den,
                            const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

void fmpz_poly_mat_solve_fflu_precomp(fmpz_poly_mat_t X,
                    const slong * perm,
                    const fmpz_poly_mat_t FFLU, const fmpz_poly_mat_t B);

#ifdef __cplusplus
}
#endif

#endif
