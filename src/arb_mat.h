/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_MAT_H
#define ARB_MAT_H

#ifdef ARB_MAT_INLINES_C
#define ARB_MAT_INLINE
#else
#define ARB_MAT_INLINE static inline
#endif

#include "fmpq_types.h"
#include "arb.h"

#ifdef __cplusplus
extern "C" {
#endif

#define arb_mat_entry(mat,i,j) ((mat)->rows[i] + (j))
#define arb_mat_nrows(mat) ((mat)->r)
#define arb_mat_ncols(mat) ((mat)->c)

ARB_MAT_INLINE arb_ptr
arb_mat_entry_ptr(arb_mat_t mat, slong i, slong j)
{
    return arb_mat_entry(mat, i, j);
}

/* Memory management */

void arb_mat_init(arb_mat_t mat, slong r, slong c);

void arb_mat_clear(arb_mat_t mat);

ARB_MAT_INLINE void
arb_mat_swap(arb_mat_t mat1, arb_mat_t mat2)
{
    FLINT_SWAP(arb_mat_struct, *mat1, *mat2);
}

ARB_MAT_INLINE void
arb_mat_swap_entrywise(arb_mat_t mat1, arb_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(mat1); i++)
        for (j = 0; j < arb_mat_ncols(mat1); j++)
            arb_swap(arb_mat_entry(mat2, i, j), arb_mat_entry(mat1, i, j));
}

/* Window matrices */

void arb_mat_window_init(arb_mat_t window, const arb_mat_t mat, slong r1, slong c1, slong r2, slong c2);

ARB_MAT_INLINE void
arb_mat_window_clear(arb_mat_t window)
{
    flint_free(window->rows);
}

/* Conversions */

void arb_mat_set(arb_mat_t dest, const arb_mat_t src);

void arb_mat_set_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src);

void arb_mat_set_round_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src, slong prec);

void arb_mat_set_fmpq_mat(arb_mat_t dest, const fmpq_mat_t src, slong prec);

/* Random generation */

void arb_mat_randtest(arb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits);

void arb_mat_randtest_cho(arb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits);

void arb_mat_randtest_spd(arb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits);

/* I/O */

#ifdef FLINT_HAVE_FILE
void arb_mat_fprintd(FILE * file, const arb_mat_t mat, slong digits);
#endif

void arb_mat_printd(const arb_mat_t mat, slong digits);

/* Comparisons */

int arb_mat_eq(const arb_mat_t mat1, const arb_mat_t mat2);

int arb_mat_ne(const arb_mat_t mat1, const arb_mat_t mat2);

int arb_mat_equal(const arb_mat_t mat1, const arb_mat_t mat2);

int arb_mat_overlaps(const arb_mat_t mat1, const arb_mat_t mat2);

int arb_mat_contains(const arb_mat_t mat1, const arb_mat_t mat2);

int arb_mat_contains_fmpq_mat(const arb_mat_t mat1, const fmpq_mat_t mat2);

int arb_mat_contains_fmpz_mat(const arb_mat_t mat1, const fmpz_mat_t mat2);

ARB_MAT_INLINE int
arb_mat_is_empty(const arb_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

ARB_MAT_INLINE int
arb_mat_is_square(const arb_mat_t mat)
{
    return (mat->r == mat->c);
}

int arb_mat_is_exact(const arb_mat_t A);

int arb_mat_is_zero(const arb_mat_t mat);
int arb_mat_is_finite(const arb_mat_t mat);
int arb_mat_is_triu(const arb_mat_t mat);
int arb_mat_is_tril(const arb_mat_t mat);

ARB_MAT_INLINE int
arb_mat_is_diag(const arb_mat_t mat)
{
    return arb_mat_is_tril(mat) && arb_mat_is_triu(mat);
}

/* Radius and interval operations */

ARB_MAT_INLINE void
arb_mat_get_mid(arb_mat_t B, const arb_mat_t A)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_get_mid_arb(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j));
}

ARB_MAT_INLINE void
arb_mat_add_error_mag(arb_mat_t mat, const mag_t err)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(mat); i++)
        for (j = 0; j < arb_mat_ncols(mat); j++)
            arb_add_error_mag(arb_mat_entry(mat, i, j), err);
}

/* Special matrices */

void arb_mat_zero(arb_mat_t mat);

void arb_mat_one(arb_mat_t mat);

void arb_mat_ones(arb_mat_t mat);

void arb_mat_indeterminate(arb_mat_t mat);

void arb_mat_hilbert(arb_mat_t mat, slong prec);

void arb_mat_pascal(arb_mat_t mat, int triangular, slong prec);

void arb_mat_stirling(arb_mat_t mat, int kind, slong prec);

void arb_mat_dct(arb_mat_t mat, int type, slong prec);

void arb_mat_transpose(arb_mat_t mat1, const arb_mat_t mat2);

/* Norms */

void arb_mat_bound_inf_norm(mag_t b, const arb_mat_t A);

void arb_mat_frobenius_norm(arb_t res, const arb_mat_t A, slong prec);

void arb_mat_bound_frobenius_norm(mag_t b, const arb_mat_t A);

/* Arithmetic */

void arb_mat_neg(arb_mat_t dest, const arb_mat_t src);

void arb_mat_add(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, slong prec);

void arb_mat_sub(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, slong prec);

void arb_mat_mul(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, slong prec);

void arb_mat_mul_classical(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec);

void arb_mat_mul_threaded(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec);

void _arb_mat_addmul_rad_mag_fast(arb_mat_t C, mag_srcptr A, mag_srcptr B, slong ar, slong ac, slong bc);

void arb_mat_mul_block(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec);

void arb_mat_mul_entrywise(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, slong prec);

void arb_mat_sqr_classical(arb_mat_t B, const arb_mat_t A, slong prec);
void arb_mat_sqr(arb_mat_t B, const arb_mat_t A, slong prec);

void arb_mat_pow_ui(arb_mat_t B, const arb_mat_t A, ulong exp, slong prec);

/* Scalar arithmetic */

ARB_MAT_INLINE void
arb_mat_scalar_mul_2exp_si(arb_mat_t B, const arb_mat_t A, slong c)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_mul_2exp_si(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c);
}

ARB_MAT_INLINE void
arb_mat_scalar_addmul_si(arb_mat_t B, const arb_mat_t A, slong c, slong prec)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_addmul_si(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_mul_si(arb_mat_t B, const arb_mat_t A, slong c, slong prec)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_mul_si(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_div_si(arb_mat_t B, const arb_mat_t A, slong c, slong prec)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_div_si(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_addmul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_addmul_fmpz(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_mul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_mul_fmpz(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_div_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_div_fmpz(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_addmul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_addmul(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_mul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_mul(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_div_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_div(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

/* Vector arithmetic */

void _arb_mat_vector_mul_row(arb_ptr res, arb_srcptr v, const arb_mat_t A, slong prec);

void _arb_mat_vector_mul_col(arb_ptr res, const arb_mat_t A, arb_srcptr v, slong prec);

void arb_mat_vector_mul_row(arb_ptr res, arb_srcptr v, const arb_mat_t A, slong prec);

void arb_mat_vector_mul_col(arb_ptr res, const arb_mat_t A, arb_srcptr v, slong prec);

/* Solving */

ARB_MAT_INLINE void
arb_mat_swap_rows(arb_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s)
    {
        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

        FLINT_SWAP(arb_ptr, mat->rows[r], mat->rows[s]);
    }
}

slong arb_mat_find_pivot_partial(const arb_mat_t mat,
                                    slong start_row, slong end_row, slong c);

void arb_mat_solve_tril_classical(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, slong prec);
void arb_mat_solve_tril_recursive(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, slong prec);
void arb_mat_solve_tril(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, slong prec);

void arb_mat_solve_triu_classical(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, slong prec);
void arb_mat_solve_triu_recursive(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, slong prec);
void arb_mat_solve_triu(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, slong prec);

int arb_mat_lu_classical(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec);
int arb_mat_lu_recursive(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec);
int arb_mat_lu(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec);

void arb_mat_solve_lu_precomp(arb_mat_t X, const slong * perm,
    const arb_mat_t A, const arb_mat_t B, slong prec);

int arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec);

int arb_mat_solve_lu(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec);

int arb_mat_solve_precond(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec);

int arb_mat_solve_preapprox(arb_mat_t X, const arb_mat_t A,
    const arb_mat_t B, const arb_mat_t R, const arb_mat_t T, slong prec);

void arb_mat_approx_mul(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec);
void arb_mat_approx_solve_triu(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, slong prec);
void arb_mat_approx_solve_tril(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, slong prec);
int arb_mat_approx_lu(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec);
void arb_mat_approx_solve_lu_precomp(arb_mat_t X, const slong * perm, const arb_mat_t A, const arb_mat_t B, slong prec);
int arb_mat_approx_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec);
int arb_mat_approx_inv(arb_mat_t X, const arb_mat_t A, slong prec);

int arb_mat_inv(arb_mat_t X, const arb_mat_t A, slong prec);

void arb_mat_det_lu(arb_t det, const arb_mat_t A, slong prec);
void arb_mat_det_precond(arb_t det, const arb_mat_t A, slong prec);
void arb_mat_det(arb_t det, const arb_mat_t A, slong prec);

int _arb_mat_cholesky_banachiewicz(arb_mat_t A, slong prec);

int arb_mat_cho(arb_mat_t L, const arb_mat_t A, slong prec);

void arb_mat_solve_cho_precomp(arb_mat_t X,
    const arb_mat_t L, const arb_mat_t B, slong prec);

void arb_mat_inv_cho_precomp(arb_mat_t X, const arb_mat_t L, slong prec);

int arb_mat_spd_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec);

int arb_mat_spd_inv(arb_mat_t X, const arb_mat_t A, slong prec);

int _arb_mat_ldl_inplace(arb_mat_t A, slong prec);

int _arb_mat_ldl_golub_and_van_loan(arb_mat_t A, slong prec);

int arb_mat_ldl(arb_mat_t L, const arb_mat_t A, slong prec);

void arb_mat_solve_ldl_precomp(arb_mat_t X,
    const arb_mat_t L, const arb_mat_t B, slong prec);

void arb_mat_inv_ldl_precomp(arb_mat_t X, const arb_mat_t L, slong prec);

/* Special functions */

void arb_mat_exp_taylor_sum(arb_mat_t S, const arb_mat_t A, slong N, slong prec);

void arb_mat_exp(arb_mat_t B, const arb_mat_t A, slong prec);

void _arb_mat_charpoly(arb_ptr poly, const arb_mat_t mat, slong prec);
void arb_mat_charpoly(arb_poly_t poly, const arb_mat_t mat, slong prec);
void _arb_mat_companion(arb_mat_t mat, arb_srcptr poly, slong prec);
void arb_mat_companion(arb_mat_t mat, const arb_poly_t poly, slong prec);

void arb_mat_trace(arb_t trace, const arb_mat_t mat, slong prec);

void _arb_mat_diag_prod(arb_t res, const arb_mat_t A, slong a, slong b, slong prec);
void arb_mat_diag_prod(arb_t res, const arb_mat_t A, slong prec);

/* Sparsity structure */

void arb_mat_entrywise_is_zero(fmpz_mat_t dest, const arb_mat_t src);

void arb_mat_entrywise_not_is_zero(fmpz_mat_t dest, const arb_mat_t src);

slong arb_mat_count_is_zero(const arb_mat_t mat);

ARB_MAT_INLINE slong
arb_mat_count_not_is_zero(const arb_mat_t mat)
{
    slong size;
    size = arb_mat_nrows(mat) * arb_mat_ncols(mat);
    return size - arb_mat_count_is_zero(mat);
}

ARB_MAT_INLINE slong
arb_mat_allocated_bytes(const arb_mat_t x)
{
    return _arb_vec_allocated_bytes(x->entries, x->r * x->c) + x->r * sizeof(arb_ptr);
}

/* LLL reduction */

int arb_mat_spd_get_fmpz_mat(fmpz_mat_t B, const arb_mat_t A, slong prec);

void arb_mat_spd_lll_reduce(fmpz_mat_t U, const arb_mat_t A, slong prec);

int arb_mat_spd_is_lll_reduced(const arb_mat_t A, slong tol_exp, slong prec);

#ifdef __cplusplus
}
#endif

#endif

