/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_MAT_H
#define ACB_MAT_H

#ifdef ACB_MAT_INLINES_C
#define ACB_MAT_INLINE
#else
#define ACB_MAT_INLINE static inline
#endif

#include "fmpq_types.h"
#include "arb_mat.h"
#include "acb.h"

#ifdef __cplusplus
extern "C" {
#endif

#define acb_mat_entry(mat,i,j) ((mat)->rows[i] + (j))
#define acb_mat_nrows(mat) ((mat)->r)
#define acb_mat_ncols(mat) ((mat)->c)

ACB_MAT_INLINE acb_ptr
acb_mat_entry_ptr(acb_mat_t mat, slong i, slong j)
{
    return acb_mat_entry(mat, i, j);
}

/* Memory management */

void acb_mat_init(acb_mat_t mat, slong r, slong c);

void acb_mat_clear(acb_mat_t mat);

ACB_MAT_INLINE void
acb_mat_swap(acb_mat_t mat1, acb_mat_t mat2)
{
    acb_mat_struct t = *mat1;
    *mat1 = *mat2;
    *mat2 = t;
}

ACB_MAT_INLINE void
acb_mat_swap_entrywise(acb_mat_t mat1, acb_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(mat1); i++)
        for (j = 0; j < acb_mat_ncols(mat1); j++)
            acb_swap(acb_mat_entry(mat2, i, j), acb_mat_entry(mat1, i, j));
}

/* Window matrices */

void acb_mat_window_init(acb_mat_t window, const acb_mat_t mat, slong r1, slong c1, slong r2, slong c2);

ARB_MAT_INLINE void
acb_mat_window_clear(acb_mat_t window)
{
    flint_free(window->rows);
}

/* Conversions */

void acb_mat_set(acb_mat_t dest, const acb_mat_t src);

void acb_mat_set_fmpz_mat(acb_mat_t dest, const fmpz_mat_t src);

void acb_mat_set_round_fmpz_mat(acb_mat_t dest, const fmpz_mat_t src, slong prec);

void acb_mat_set_fmpq_mat(acb_mat_t dest, const fmpq_mat_t src, slong prec);

void acb_mat_set_arb_mat(acb_mat_t dest, const arb_mat_t src);

void acb_mat_set_round_arb_mat(acb_mat_t dest, const arb_mat_t src, slong prec);

void acb_mat_get_real(arb_mat_t re, const acb_mat_t mat);

void acb_mat_get_imag(arb_mat_t im, const acb_mat_t mat);

void acb_mat_set_real_imag(acb_mat_t mat, const arb_mat_t re, const arb_mat_t im);

/* Random generation */

void acb_mat_randtest(acb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits);

void acb_mat_randtest_eig(acb_mat_t A, flint_rand_t state, acb_srcptr E, slong prec);

/* I/O */

#ifdef FLINT_HAVE_FILE
void acb_mat_fprintd(FILE * file, const acb_mat_t mat, slong digits);
#endif

void acb_mat_printd(const acb_mat_t mat, slong digits);

/* Comparisons */

int acb_mat_eq(const acb_mat_t mat1, const acb_mat_t mat2);

int acb_mat_ne(const acb_mat_t mat1, const acb_mat_t mat2);

int acb_mat_equal(const acb_mat_t mat1, const acb_mat_t mat2);

int acb_mat_overlaps(const acb_mat_t mat1, const acb_mat_t mat2);

int acb_mat_contains(const acb_mat_t mat1, const acb_mat_t mat2);

int acb_mat_contains_fmpq_mat(const acb_mat_t mat1, const fmpq_mat_t mat2);

int acb_mat_contains_fmpz_mat(const acb_mat_t mat1, const fmpz_mat_t mat2);

int acb_mat_is_real(const acb_mat_t mat);

ACB_MAT_INLINE int
acb_mat_is_empty(const acb_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

ACB_MAT_INLINE int
acb_mat_is_square(const acb_mat_t mat)
{
    return (mat->r == mat->c);
}

int acb_mat_is_exact(const acb_mat_t mat);

int acb_mat_is_zero(const acb_mat_t mat);
int acb_mat_is_finite(const acb_mat_t mat);
int acb_mat_is_triu(const acb_mat_t mat);
int acb_mat_is_tril(const acb_mat_t mat);

ACB_MAT_INLINE int
acb_mat_is_diag(const acb_mat_t mat)
{
    return acb_mat_is_tril(mat) && acb_mat_is_triu(mat);
}

/* Radius and interval operations */

ACB_MAT_INLINE void
acb_mat_get_mid(acb_mat_t B, const acb_mat_t A)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
    {
        for (j = 0; j < acb_mat_ncols(A); j++)
        {
            arb_get_mid_arb(acb_realref(acb_mat_entry(B, i, j)), acb_realref(acb_mat_entry(A, i, j)));
            arb_get_mid_arb(acb_imagref(acb_mat_entry(B, i, j)), acb_imagref(acb_mat_entry(A, i, j)));
        }
    }
}

ACB_MAT_INLINE void
acb_mat_add_error_mag(acb_mat_t mat, const mag_t err)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(mat); i++)
        for (j = 0; j < acb_mat_ncols(mat); j++)
            acb_add_error_mag(acb_mat_entry(mat, i, j), err);
}

/* Special matrices */

void acb_mat_zero(acb_mat_t mat);

void acb_mat_one(acb_mat_t mat);

void acb_mat_ones(acb_mat_t mat);

void acb_mat_onei(acb_mat_t mat);

void acb_mat_indeterminate(acb_mat_t mat);

void acb_mat_dft(acb_mat_t res, int kind, slong prec);

void acb_mat_transpose(acb_mat_t mat1, const acb_mat_t mat2);

void acb_mat_conjugate(acb_mat_t mat1, const acb_mat_t mat2);

ACB_MAT_INLINE void
acb_mat_conjugate_transpose(acb_mat_t mat1, const acb_mat_t mat2)
{
    acb_mat_transpose(mat1, mat2);
    acb_mat_conjugate(mat1, mat1);
}

/* Norms */

void acb_mat_bound_inf_norm(mag_t b, const acb_mat_t A);

void acb_mat_frobenius_norm(arb_t res, const acb_mat_t A, slong prec);

void acb_mat_bound_frobenius_norm(mag_t b, const acb_mat_t A);

/* Arithmetic */

void acb_mat_neg(acb_mat_t dest, const acb_mat_t src);

void acb_mat_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec);

void acb_mat_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec);

void acb_mat_mul_classical(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec);
void acb_mat_mul_threaded(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec);
void acb_mat_mul_reorder(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec);
void acb_mat_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec);

void acb_mat_mul_entrywise(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec);

void acb_mat_sqr_classical(acb_mat_t res, const acb_mat_t mat, slong prec);
void acb_mat_sqr(acb_mat_t res, const acb_mat_t mat, slong prec);

void acb_mat_pow_ui(acb_mat_t B, const acb_mat_t A, ulong exp, slong prec);

/* Scalar arithmetic */

ACB_MAT_INLINE void
acb_mat_scalar_mul_2exp_si(acb_mat_t B, const acb_mat_t A, slong c)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul_2exp_si(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c);
}

ACB_MAT_INLINE void
acb_mat_scalar_addmul_si(acb_mat_t B, const acb_mat_t A, slong c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_addmul_si(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_mul_si(acb_mat_t B, const acb_mat_t A, slong c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul_si(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_div_si(acb_mat_t B, const acb_mat_t A, slong c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_div_si(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_addmul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_addmul_fmpz(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_mul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul_fmpz(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_div_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_div_fmpz(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_addmul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_addmul(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_mul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_div_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_div(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_addmul_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_addmul_arb(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_mul_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul_arb(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_div_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, slong prec)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_div_arb(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

/* Vector arithmetic */

void _acb_mat_vector_mul_row(acb_ptr res, acb_srcptr v, const acb_mat_t A, slong prec);

void _acb_mat_vector_mul_col(acb_ptr res, const acb_mat_t A, acb_srcptr v, slong prec);

void acb_mat_vector_mul_row(acb_ptr res, acb_srcptr v, const acb_mat_t A, slong prec);

void acb_mat_vector_mul_col(acb_ptr res, const acb_mat_t A, acb_srcptr v, slong prec);

/* Solving */

ACB_MAT_INLINE void
acb_mat_swap_rows(acb_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s)
    {
        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

        FLINT_SWAP(acb_ptr, mat->rows[r], mat->rows[s]);
    }
}

slong acb_mat_find_pivot_partial(const acb_mat_t mat,
                                    slong start_row, slong end_row, slong c);

void acb_mat_solve_tril_classical(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec);
void acb_mat_solve_tril_recursive(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec);
void acb_mat_solve_tril(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec);

void acb_mat_solve_triu_classical(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec);
void acb_mat_solve_triu_recursive(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec);
void acb_mat_solve_triu(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec);

int acb_mat_lu_classical(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec);
int acb_mat_lu_recursive(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec);
int acb_mat_lu(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec);

void acb_mat_solve_lu_precomp(acb_mat_t X, const slong * perm,
    const acb_mat_t A, const acb_mat_t B, slong prec);

int acb_mat_solve_lu(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec);

int acb_mat_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec);

int acb_mat_solve_precond(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec);

void acb_mat_approx_mul(acb_mat_t C, const acb_mat_t A, const acb_mat_t B, slong prec);
void acb_mat_approx_solve_triu(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec);
void acb_mat_approx_solve_tril(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec);
int acb_mat_approx_lu(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec);
void acb_mat_approx_solve_lu_precomp(acb_mat_t X, const slong * perm, const acb_mat_t A, const acb_mat_t B, slong prec);
int acb_mat_approx_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec);
int acb_mat_approx_inv(acb_mat_t X, const acb_mat_t A, slong prec);

int acb_mat_inv(acb_mat_t X, const acb_mat_t A, slong prec);

void acb_mat_det_lu(acb_t det, const acb_mat_t A, slong prec);
void acb_mat_det_precond(acb_t det, const acb_mat_t A, slong prec);
void acb_mat_det(acb_t det, const acb_mat_t A, slong prec);

/* Eigenvalues and eigenvectors */

int acb_mat_approx_eig_qr(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, const mag_t tol, slong maxiter, slong prec);

void acb_mat_eig_global_enclosure(mag_t eps, const acb_mat_t A, acb_srcptr E, const acb_mat_t R, slong prec);

void acb_mat_eig_enclosure_rump(acb_t lambda, acb_mat_t J, acb_mat_t X, const acb_mat_t A,
    const acb_t lambda_approx, const acb_mat_t X_approx, slong prec);

int acb_mat_eig_simple_rump(acb_ptr E, acb_mat_t L, acb_mat_t R,
    const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec);
int acb_mat_eig_simple_vdhoeven_mourrain(acb_ptr E, acb_mat_t L, acb_mat_t R,
    const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec);
int acb_mat_eig_simple(acb_ptr E, acb_mat_t L, acb_mat_t R,
    const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec);

int acb_mat_eig_multiple_rump(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec);
int acb_mat_eig_multiple(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec);

/* Special functions */

void acb_mat_exp_taylor_sum(acb_mat_t S, const acb_mat_t A, slong N, slong prec);

void acb_mat_exp(acb_mat_t B, const acb_mat_t A, slong prec);

void _acb_mat_charpoly(acb_ptr poly, const acb_mat_t mat, slong prec);
void acb_mat_charpoly(acb_poly_t poly, const acb_mat_t mat, slong prec);
void _acb_mat_companion(acb_mat_t mat, acb_srcptr poly, slong prec);
void acb_mat_companion(acb_mat_t mat, const acb_poly_t poly, slong prec);

void acb_mat_trace(acb_t trace, const acb_mat_t mat, slong prec);

void _acb_mat_diag_prod(acb_t res, const acb_mat_t A, slong a, slong b, slong prec);
void acb_mat_diag_prod(acb_t res, const acb_mat_t A, slong prec);

ACB_MAT_INLINE slong
acb_mat_allocated_bytes(const acb_mat_t x)
{
    return _acb_vec_allocated_bytes(x->entries, x->r * x->c) + x->r * sizeof(acb_ptr);
}

#ifdef __cplusplus
}
#endif

#endif

