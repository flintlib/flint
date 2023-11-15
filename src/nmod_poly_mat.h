/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_H
#define NMOD_POLY_MAT_H

#ifdef NMOD_POLY_MAT_INLINES_C
#define NMOD_POLY_MAT_INLINE
#else
#define NMOD_POLY_MAT_INLINE static inline
#endif

#include "nmod_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

NMOD_POLY_MAT_INLINE
nmod_poly_struct * nmod_poly_mat_entry(const nmod_poly_mat_t mat, slong i, slong j)
{
    return mat->rows[i] + j;
}

NMOD_POLY_MAT_INLINE slong
nmod_poly_mat_nrows(const nmod_poly_mat_t mat)
{
    return mat->r;
}

NMOD_POLY_MAT_INLINE slong
nmod_poly_mat_ncols(const nmod_poly_mat_t mat)
{
    return mat->c;
}

/* Memory management *********************************************************/

void nmod_poly_mat_init(nmod_poly_mat_t mat, slong rows, slong cols, mp_limb_t n);

void nmod_poly_mat_init_set(nmod_poly_mat_t mat, const nmod_poly_mat_t src);

void nmod_poly_mat_swap(nmod_poly_mat_t mat1, nmod_poly_mat_t mat2);

NMOD_POLY_MAT_INLINE void
nmod_poly_mat_swap_entrywise(nmod_poly_mat_t mat1, nmod_poly_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < nmod_poly_mat_nrows(mat1); i++)
        for (j = 0; j < nmod_poly_mat_ncols(mat1); j++)
            nmod_poly_swap(nmod_poly_mat_entry(mat2, i, j), nmod_poly_mat_entry(mat1, i, j));
}

void nmod_poly_mat_set(nmod_poly_mat_t mat1, const nmod_poly_mat_t mat2);

void nmod_poly_mat_set_nmod_mat(nmod_poly_mat_t pmat, const nmod_mat_t cmat);

void nmod_poly_mat_clear(nmod_poly_mat_t mat);

/* Truncate, shift *********************************************************/

void nmod_poly_mat_set_trunc(nmod_poly_mat_t res,
                             const nmod_poly_mat_t pmat,
                             long len);

NMOD_POLY_MAT_INLINE
void nmod_poly_mat_truncate(nmod_poly_mat_t pmat, long len)
{
    nmod_poly_mat_set_trunc(pmat, pmat, len);
}

void nmod_poly_mat_shift_left(nmod_poly_mat_t res,
                              const nmod_poly_mat_t pmat,
                              slong k);

void nmod_poly_mat_shift_right(nmod_poly_mat_t res,
                               const nmod_poly_mat_t pmat,
                               slong k);

/* Basic properties **********************************************************/

NMOD_POLY_MAT_INLINE mp_limb_t
nmod_poly_mat_modulus(const nmod_poly_mat_t mat)
{
    return mat->modulus;
}

void nmod_poly_mat_get_coeff_mat(nmod_mat_t coeff,
                                 const nmod_poly_mat_t pmat,
                                 slong deg);

void nmod_poly_mat_set_coeff_mat(nmod_poly_mat_t pmat,
                                 const nmod_mat_t coeff,
                                 slong deg);


/* Comparison ****************************************************************/

int nmod_poly_mat_equal(const nmod_poly_mat_t mat1,
                        const nmod_poly_mat_t mat2);

int nmod_poly_mat_equal_nmod_mat(const nmod_poly_mat_t pmat,
                                const nmod_mat_t cmat);

int nmod_poly_mat_is_zero(const nmod_poly_mat_t mat);

int nmod_poly_mat_is_one(const nmod_poly_mat_t mat);

NMOD_POLY_MAT_INLINE int
nmod_poly_mat_is_empty(const nmod_poly_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

NMOD_POLY_MAT_INLINE int
nmod_poly_mat_is_square(const nmod_poly_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Standard matrices *********************************************************/

void nmod_poly_mat_zero(nmod_poly_mat_t mat);

void nmod_poly_mat_one(nmod_poly_mat_t mat);

/* Random matrices ***********************************************************/

void nmod_poly_mat_randtest(nmod_poly_mat_t mat, flint_rand_t state,
                                slong len);

void nmod_poly_mat_randtest_sparse(nmod_poly_mat_t A, flint_rand_t state,
                        slong len, float density);

/* Windows and concatenation */

void nmod_poly_mat_window_init(nmod_poly_mat_t window, const nmod_poly_mat_t mat, slong r1,
    slong c1, slong r2, slong c2);

void nmod_poly_mat_window_clear(nmod_poly_mat_t window);

void nmod_poly_mat_concat_horizontal(nmod_poly_mat_t res,
                           const nmod_poly_mat_t mat1,  const nmod_poly_mat_t mat2);

void nmod_poly_mat_concat_vertical(nmod_poly_mat_t res,
                           const nmod_poly_mat_t mat1,  const nmod_poly_mat_t mat2);

/* Input and output **********************************************************/

void nmod_poly_mat_print(const nmod_poly_mat_t mat, const char * x);

/* Norms *********************************************************************/

slong nmod_poly_mat_max_length(const nmod_poly_mat_t A);

NMOD_POLY_MAT_INLINE
slong nmod_poly_mat_degree(const nmod_poly_mat_t pmat)
{
    return nmod_poly_mat_max_length(pmat)-1;
}


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

slong nmod_poly_mat_find_pivot_any(const nmod_poly_mat_t mat,
                                    slong start_row, slong end_row, slong c);

slong nmod_poly_mat_find_pivot_partial(const nmod_poly_mat_t mat,
                                    slong start_row, slong end_row, slong c);

slong nmod_poly_mat_fflu(nmod_poly_mat_t B, nmod_poly_t den, slong * perm,
                            const nmod_poly_mat_t A, int rank_check);

slong nmod_poly_mat_rref(nmod_poly_mat_t B, nmod_poly_t den,
                            const nmod_poly_mat_t A);

/* Trace *********************************************************************/

void nmod_poly_mat_trace(nmod_poly_t trace, const nmod_poly_mat_t mat);

/* Determinant and rank ******************************************************/

void nmod_poly_mat_det(nmod_poly_t det, const nmod_poly_mat_t A);

void nmod_poly_mat_det_fflu(nmod_poly_t det, const nmod_poly_mat_t A);

void nmod_poly_mat_det_interpolate(nmod_poly_t det, const nmod_poly_mat_t A);

slong nmod_poly_mat_rank(const nmod_poly_mat_t A);

/* Inverse *******************************************************************/

int nmod_poly_mat_inv(nmod_poly_mat_t Ainv, nmod_poly_t den,
    const nmod_poly_mat_t A);

/* Nullspace *****************************************************************/

slong nmod_poly_mat_nullspace(nmod_poly_mat_t res, const nmod_poly_mat_t mat);

/* Solving *******************************************************************/

int nmod_poly_mat_solve(nmod_poly_mat_t X, nmod_poly_t den,
                    const nmod_poly_mat_t A, const nmod_poly_mat_t B);

int nmod_poly_mat_solve_fflu(nmod_poly_mat_t X, nmod_poly_t den,
                            const nmod_poly_mat_t A, const nmod_poly_mat_t B);

void nmod_poly_mat_solve_fflu_precomp(nmod_poly_mat_t X,
                    const slong * perm,
                    const nmod_poly_mat_t FFLU, const nmod_poly_mat_t B);

#ifdef __cplusplus
}
#endif

#endif
