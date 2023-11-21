/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_MAT_H
#define CA_MAT_H

#ifdef CA_MAT_INLINES_C
#define CA_MAT_INLINE
#else
#define CA_MAT_INLINE static inline
#endif

#include "ca_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Matrix object */

typedef struct
{
    ca_ptr entries;
    slong r;
    slong c;
    ca_ptr * rows;
}
ca_mat_struct;

typedef ca_mat_struct ca_mat_t[1];

#define ca_mat_entry(mat,i,j) ((mat)->rows[i] + (j))
#define ca_mat_nrows(mat) ((mat)->r)
#define ca_mat_ncols(mat) ((mat)->c)

CA_MAT_INLINE ca_ptr
ca_mat_entry_ptr(ca_mat_t mat, slong i, slong j)
{
    return ca_mat_entry(mat, i, j);
}

/* Memory management */

void ca_mat_init(ca_mat_t mat, slong r, slong c, ca_ctx_t ctx);

void ca_mat_clear(ca_mat_t mat, ca_ctx_t ctx);

CA_MAT_INLINE void
ca_mat_swap(ca_mat_t mat1, ca_mat_t mat2, ca_ctx_t ctx)
{
    FLINT_SWAP(ca_mat_struct, *mat1, *mat2);
}

/* Window matrices */

void ca_mat_window_init(ca_mat_t window, const ca_mat_t mat, slong r1, slong c1, slong r2, slong c2, ca_ctx_t ctx);

CA_MAT_INLINE void
ca_mat_window_clear(ca_mat_t window, ca_ctx_t ctx)
{
    flint_free(window->rows);
}

/* Shape */

CA_MAT_INLINE int
ca_mat_is_empty(const ca_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

CA_MAT_INLINE int
ca_mat_is_square(const ca_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Conversions */

void ca_mat_set(ca_mat_t dest, const ca_mat_t src, ca_ctx_t ctx);
void ca_mat_set_fmpz_mat(ca_mat_t dest, const fmpz_mat_t src, ca_ctx_t ctx);
void ca_mat_set_fmpq_mat(ca_mat_t dest, const fmpq_mat_t src, ca_ctx_t ctx);
void ca_mat_set_ca(ca_mat_t y, const ca_t x, ca_ctx_t ctx);

void ca_mat_transfer(ca_mat_t res, ca_ctx_t res_ctx, const ca_mat_t src, ca_ctx_t src_ctx);

/* Random generation */

void ca_mat_randtest(ca_mat_t mat, flint_rand_t state, slong len, slong bits, ca_ctx_t ctx);
void ca_mat_randtest_rational(ca_mat_t mat, flint_rand_t state, slong bits, ca_ctx_t ctx);
void ca_mat_randops(ca_mat_t mat, flint_rand_t state, slong count, ca_ctx_t ctx);

/* I/O */

void ca_mat_print(const ca_mat_t mat, ca_ctx_t ctx);
void ca_mat_printn(const ca_mat_t mat, slong digits, ca_ctx_t ctx);

/* Special matrices */

void ca_mat_zero(ca_mat_t mat, ca_ctx_t ctx);
void ca_mat_one(ca_mat_t mat, ca_ctx_t ctx);
void ca_mat_ones(ca_mat_t mat, ca_ctx_t ctx);
void ca_mat_pascal(ca_mat_t mat, int triangular, ca_ctx_t ctx);
void ca_mat_stirling(ca_mat_t mat, int kind, ca_ctx_t ctx);
void ca_mat_hilbert(ca_mat_t mat, ca_ctx_t ctx);
void ca_mat_dft(ca_mat_t res, int type, ca_ctx_t ctx);

/* Comparisons and properties */

truth_t ca_mat_check_equal(const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx);
truth_t ca_mat_check_is_zero(const ca_mat_t A, ca_ctx_t ctx);
truth_t ca_mat_check_is_one(const ca_mat_t A, ca_ctx_t ctx);

/* Conjugate and transpose */

void ca_mat_transpose(ca_mat_t B, const ca_mat_t A, ca_ctx_t ctx);
void ca_mat_conj(ca_mat_t B, const ca_mat_t A, ca_ctx_t ctx);
void ca_mat_conj_transpose(ca_mat_t mat1, const ca_mat_t mat2, ca_ctx_t ctx);

/* Arithmetic */

void ca_mat_add_ca(ca_mat_t y, const ca_mat_t a, const ca_t x, ca_ctx_t ctx);
void ca_mat_sub_ca(ca_mat_t y, const ca_mat_t a, const ca_t x, ca_ctx_t ctx);
void ca_mat_addmul_ca(ca_mat_t y, const ca_mat_t a, const ca_t x, ca_ctx_t ctx);
void ca_mat_submul_ca(ca_mat_t y, const ca_mat_t a, const ca_t x, ca_ctx_t ctx);

void ca_mat_neg(ca_mat_t dest, const ca_mat_t src, ca_ctx_t ctx);
void ca_mat_add(ca_mat_t res, const ca_mat_t mat1, const ca_mat_t mat2, ca_ctx_t ctx);
void ca_mat_sub(ca_mat_t res, const ca_mat_t mat1, const ca_mat_t mat2, ca_ctx_t ctx);

void ca_mat_mul(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx);

void ca_mat_mul_classical(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx);
void ca_mat_mul_same_nf(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_field_t K, ca_ctx_t ctx);

CA_MAT_INLINE void
ca_mat_mul_si(ca_mat_t B, const ca_mat_t A, slong c, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            ca_mul_si(ca_mat_entry(B, i, j), ca_mat_entry(A, i, j), c, ctx);
}

CA_MAT_INLINE void
ca_mat_mul_fmpz(ca_mat_t B, const ca_mat_t A, const fmpz_t c, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            ca_mul_fmpz(ca_mat_entry(B, i, j), ca_mat_entry(A, i, j), c, ctx);
}

CA_MAT_INLINE void
ca_mat_mul_fmpq(ca_mat_t B, const ca_mat_t A, const fmpq_t c, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            ca_mul_fmpq(ca_mat_entry(B, i, j), ca_mat_entry(A, i, j), c, ctx);
}

CA_MAT_INLINE void
ca_mat_mul_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            ca_mul(ca_mat_entry(B, i, j), ca_mat_entry(A, i, j), c, ctx);
}

CA_MAT_INLINE void
ca_mat_div_si(ca_mat_t B, const ca_mat_t A, slong c, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            ca_div_si(ca_mat_entry(B, i, j), ca_mat_entry(A, i, j), c, ctx);
}

CA_MAT_INLINE void
ca_mat_div_fmpz(ca_mat_t B, const ca_mat_t A, const fmpz_t c, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            ca_div_fmpz(ca_mat_entry(B, i, j), ca_mat_entry(A, i, j), c, ctx);
}

CA_MAT_INLINE void
ca_mat_div_fmpq(ca_mat_t B, const ca_mat_t A, const fmpq_t c, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            ca_div_fmpq(ca_mat_entry(B, i, j), ca_mat_entry(A, i, j), c, ctx);
}

CA_MAT_INLINE void
ca_mat_div_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            ca_div(ca_mat_entry(B, i, j), ca_mat_entry(A, i, j), c, ctx);
}

CA_MAT_INLINE void
ca_mat_sqr(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_mat_mul(res, A, A, ctx);
}

void ca_mat_pow_ui_binexp(ca_mat_t B, const ca_mat_t A, ulong exp, ca_ctx_t ctx);

/* Polynomial evaluation */

void _ca_mat_ca_poly_evaluate(ca_mat_t y, ca_srcptr poly, slong len, const ca_mat_t x, ca_ctx_t ctx);
void ca_mat_ca_poly_evaluate(ca_mat_t res, const ca_poly_t f, const ca_mat_t a, ca_ctx_t ctx);

/* Trace */

void ca_mat_trace(ca_t trace, const ca_mat_t mat, ca_ctx_t ctx);

/* Gaussian elimination, solving and inverse */

truth_t ca_mat_find_pivot(slong * pivot_row, ca_mat_t mat, slong start_row, slong end_row, slong column, ca_ctx_t ctx);

CA_MAT_INLINE void
_ca_mat_swap_rows(ca_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s)
    {
        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

        FLINT_SWAP(ca_ptr, mat->rows[r], mat->rows[s]);
    }
}

int ca_mat_lu_classical(slong * rank, slong * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx);
int ca_mat_lu_recursive(slong * rank, slong * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx);
int ca_mat_lu(slong * rank, slong * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx);

int ca_mat_fflu(slong * rank, slong * P, ca_mat_t LU, ca_t den, const ca_mat_t A, int rank_check, ca_ctx_t ctx);

int ca_mat_rref_fflu(slong * rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx);
int ca_mat_rref_lu(slong * rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx);
int ca_mat_rref(slong * rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx);

truth_t ca_mat_nonsingular_lu(slong * P, ca_mat_t LU, const ca_mat_t A, ca_ctx_t ctx);
truth_t ca_mat_nonsingular_fflu(slong * P, ca_mat_t LU, ca_t den, const ca_mat_t A, ca_ctx_t ctx);

truth_t ca_mat_nonsingular_solve_adjugate(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx);
truth_t ca_mat_nonsingular_solve_fflu(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx);
truth_t ca_mat_nonsingular_solve_lu(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx);
truth_t ca_mat_nonsingular_solve(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx);

truth_t ca_mat_inv(ca_mat_t X, const ca_mat_t A, ca_ctx_t ctx);

void ca_mat_solve_tril_classical(ca_mat_t X, const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx);
void ca_mat_solve_tril_recursive(ca_mat_t X, const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx);
void ca_mat_solve_tril(ca_mat_t X, const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx);

void ca_mat_solve_triu_classical(ca_mat_t X, const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx);
void ca_mat_solve_triu_recursive(ca_mat_t X, const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx);
void ca_mat_solve_triu(ca_mat_t X, const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx);

void ca_mat_solve_lu_precomp(ca_mat_t X, const slong * perm, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx);
void ca_mat_solve_fflu_precomp(ca_mat_t X, const slong * perm, const ca_mat_t A, const ca_t den, const ca_mat_t B, ca_ctx_t ctx);

/* Rank and kernel */

int ca_mat_rank(slong * rank, const ca_mat_t A, ca_ctx_t ctx);

int ca_mat_right_kernel(ca_mat_t X, const ca_mat_t A, ca_ctx_t ctx);

/* Determinant */

void ca_mat_det_berkowitz(ca_t det, const ca_mat_t A, ca_ctx_t ctx);
int ca_mat_det_lu(ca_t det, const ca_mat_t A, ca_ctx_t ctx);
int ca_mat_det_bareiss(ca_t det, const ca_mat_t A, ca_ctx_t ctx);
void ca_mat_det_cofactor(ca_t det, const ca_mat_t A, ca_ctx_t ctx);
void ca_mat_det(ca_t det, const ca_mat_t A, ca_ctx_t ctx);

void ca_mat_adjugate_cofactor(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx);
void ca_mat_adjugate_charpoly(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx);
void ca_mat_adjugate(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx);

/* Characteristic polynomial */

void _ca_mat_charpoly_berkowitz(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx);
void ca_mat_charpoly_berkowitz(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx);

int _ca_mat_charpoly_danilevsky(ca_ptr p, const ca_mat_t A, ca_ctx_t ctx);
int ca_mat_charpoly_danilevsky(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx);

void _ca_mat_charpoly(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx);
void ca_mat_charpoly(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx);



int ca_mat_companion(ca_mat_t A, const ca_poly_t poly, ca_ctx_t ctx);

/* Eigenvalues and eigenvectors */

int ca_mat_eigenvalues(ca_vec_t lambda, ulong * exp, const ca_mat_t mat, ca_ctx_t ctx);

/* Diagonalization */

truth_t ca_mat_diagonalization(ca_mat_t D, ca_mat_t P, const ca_mat_t A, ca_ctx_t ctx);

void ca_mat_set_jordan_blocks(ca_mat_t mat, const ca_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, ca_ctx_t ctx);
int ca_mat_jordan_blocks(ca_vec_t lambda, slong * num_blocks, slong * block_lambda, slong * block_size, const ca_mat_t A, ca_ctx_t ctx);
int ca_mat_jordan_transformation(ca_mat_t mat, const ca_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, const ca_mat_t A, ca_ctx_t ctx);
int ca_mat_jordan_form(ca_mat_t J, ca_mat_t P, const ca_mat_t A, ca_ctx_t ctx);

/* Matrix functions */

int ca_mat_exp(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx);
truth_t ca_mat_log(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx);

/* Internal representation */

/* todo: document, make consistent */
ca_field_ptr _ca_mat_same_field(const ca_mat_t A, ca_ctx_t ctx);
ca_field_ptr _ca_mat_same_field2(const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
