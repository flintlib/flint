/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_MAT_H
#define GR_MAT_H

#ifdef GR_MAT_INLINES_C
#define GR_MAT_INLINE
#else
#define GR_MAT_INLINE static inline
#endif

#include "fmpq_types.h"
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define GR_MAT_ENTRY(mat,i,j,sz) GR_ENTRY((mat)->rows[i], j, sz)
#define gr_mat_nrows(mat, ctx) ((mat)->r)
#define gr_mat_ncols(mat, ctx) ((mat)->c)

GR_MAT_INLINE gr_ptr gr_mat_entry_ptr(gr_mat_t mat, slong i, slong j, gr_ctx_t ctx)
{
    return GR_MAT_ENTRY(mat, i, j, ctx->sizeof_elem);
}

GR_MAT_INLINE gr_srcptr gr_mat_entry_srcptr(const gr_mat_t mat, slong i, slong j, gr_ctx_t ctx)
{
    return GR_MAT_ENTRY(mat, i, j, ctx->sizeof_elem);
}

/* Generics */
typedef int ((*gr_method_mat_unary_op_get_scalar)(gr_ptr, const gr_mat_t, gr_ctx_ptr));
typedef int ((*gr_method_mat_unary_op)(gr_mat_t, const gr_mat_t, gr_ctx_ptr));
typedef int ((*gr_method_mat_binary_op)(gr_mat_t, const gr_mat_t, const gr_mat_t, gr_ctx_ptr));
typedef int ((*gr_method_mat_binary_op_with_flag)(gr_mat_t, const gr_mat_t, const gr_mat_t, int, gr_ctx_ptr));
typedef int ((*gr_method_mat_pivot_op)(slong *, gr_mat_t, slong, slong, slong, gr_ctx_ptr));
typedef int ((*gr_method_mat_diagonalization_op)(gr_vec_t, gr_mat_t, gr_mat_t, const gr_mat_t, int, gr_ctx_ptr));
typedef int ((*gr_method_mat_lu_op)(slong *, slong *, gr_mat_t, const gr_mat_t, int, gr_ctx_ptr));

#define GR_MAT_UNARY_OP_GET_SCALAR(ctx, NAME) (((gr_method_mat_unary_op_get_scalar *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_MAT_UNARY_OP(ctx, NAME) (((gr_method_mat_unary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_MAT_BINARY_OP(ctx, NAME) (((gr_method_mat_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_MAT_BINARY_OP_WITH_FLAG(ctx, NAME) (((gr_method_mat_binary_op_with_flag *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_MAT_PIVOT_OP(ctx, NAME) (((gr_method_mat_pivot_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_MAT_DIAGONALIZATION_OP(ctx, NAME) (((gr_method_mat_diagonalization_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_MAT_LU_OP(ctx, NAME) (((gr_method_mat_lu_op *) ctx->methods)[GR_METHOD_ ## NAME])

void gr_mat_init(gr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_init_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
void gr_mat_clear(gr_mat_t mat, gr_ctx_t ctx);

GR_MAT_INLINE void
gr_mat_swap(gr_mat_t mat1, gr_mat_t mat2, gr_ctx_t FLINT_UNUSED(ctx))
{
    FLINT_SWAP(gr_mat_struct, *mat1, *mat2);
}

WARN_UNUSED_RESULT int gr_mat_swap_rows(gr_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_invert_rows(gr_mat_t mat, slong * perm, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_swap_cols(gr_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_invert_cols(gr_mat_t mat, slong * perm, gr_ctx_t ctx);

void gr_mat_window_init(gr_mat_t window, const gr_mat_t mat, slong r1, slong c1, slong r2, slong c2, gr_ctx_t ctx);

GR_MAT_INLINE void
gr_mat_window_clear(gr_mat_t window, gr_ctx_t FLINT_UNUSED(ctx))
{
    flint_free(window->rows);
}

WARN_UNUSED_RESULT int gr_mat_concat_horizontal(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_concat_vertical(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);

int gr_mat_write(gr_stream_t out, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_print(const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_randtest(gr_mat_t mat, flint_rand_t state, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_randops(gr_mat_t mat, flint_rand_t state, slong count, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_randpermdiag(int * parity, gr_mat_t mat, flint_rand_t state, gr_ptr diag, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_randrank(gr_mat_t mat, flint_rand_t state, slong rank, gr_ctx_t ctx);

GR_MAT_INLINE truth_t
gr_mat_is_empty(const gr_mat_t mat, gr_ctx_t FLINT_UNUSED(ctx))
{
    return ((mat->r == 0) || (mat->c == 0)) ? T_TRUE : T_FALSE;
}

GR_MAT_INLINE truth_t
gr_mat_is_square(const gr_mat_t mat, gr_ctx_t FLINT_UNUSED(ctx))
{
    return (mat->r == mat->c) ? T_TRUE : T_FALSE;
}

truth_t gr_mat_equal(const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
truth_t gr_mat_is_zero(const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_one(const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_zero(gr_mat_t res, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_one(gr_mat_t res, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_set_scalar(gr_mat_t res, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_set_ui(gr_mat_t res, ulong v, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_set_si(gr_mat_t res, slong v, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_set_fmpz(gr_mat_t res, const fmpz_t v, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_set_fmpq(gr_mat_t res, const fmpq_t v, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_set_fmpz_mat(gr_mat_t res, const fmpz_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_set_fmpq_mat(gr_mat_t res, const fmpq_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_set_gr_mat_other(gr_mat_t res, const gr_mat_t mat, gr_ctx_t mat_ctx, gr_ctx_t ctx);

/* fixme: needed for method typedefs */
#ifdef GR_H
WARN_UNUSED_RESULT int gr_mat_entrywise_unary_op(gr_mat_t res, gr_method_unary_op f, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_entrywise_binary_op(gr_mat_t res, gr_method_binary_op f, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_entrywise_binary_op_scalar(gr_mat_t res, gr_method_binary_op f, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx);
truth_t gr_mat_entrywise_unary_predicate_all(gr_method_unary_predicate f, const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_entrywise_unary_predicate_any(gr_method_unary_predicate f, const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_entrywise_binary_predicate_all(gr_method_binary_predicate f, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
#endif

WARN_UNUSED_RESULT int gr_mat_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_swap_entrywise(gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);

/* todo: test, wrap; div; more conversions */
WARN_UNUSED_RESULT int gr_mat_add_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_sub_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_mul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_addmul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_submul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_div_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_mul_classical(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_mul_strassen(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_mul_generic(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

/* todo */
GR_MAT_INLINE WARN_UNUSED_RESULT int
gr_mat_sqr(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_mul(res, mat, mat, ctx);
}

WARN_UNUSED_RESULT int _gr_mat_gr_poly_evaluate(gr_mat_t y, gr_srcptr poly, slong len, const gr_mat_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_gr_poly_evaluate(gr_mat_t res, const gr_poly_t f, const gr_mat_t a, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_find_nonzero_pivot_large_abs(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_find_nonzero_pivot_generic(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_find_nonzero_pivot(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_lu_recursive(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_lu_classical(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_lu_generic(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_fflu(slong * res_rank, slong * P, gr_mat_t LU, gr_ptr den, const gr_mat_t A, int rank_check, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_fflu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_lu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_fflu_precomp(gr_mat_t X, const slong * perm, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_lu_precomp(gr_mat_t X, const slong * perm, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_den_fflu(gr_mat_t X, gr_ptr den, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_den(gr_mat_t X, gr_ptr den, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_solve_field(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_det_berkowitz(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_det_fflu(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_det_lu(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_det_cofactor(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_det_generic_field(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_det_generic_integral_domain(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_det_generic(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_det(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_inv(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_adjugate_charpoly(gr_mat_t adj, gr_ptr det, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_adjugate_cofactor(gr_mat_t adj, gr_ptr det, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_adjugate(gr_mat_t adj, gr_ptr det, const gr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_rank_lu(slong * rank, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_rank_fflu(slong * rank, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_rank(slong * rank, const gr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_rref_lu(slong * res_rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_rref_fflu(slong * res_rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_rref(slong * res_rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_rref_den_fflu(slong * res_rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_rref_den(slong * res_rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_nullspace(gr_mat_t X, const gr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_ones(gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_pascal(gr_mat_t mat, int triangular, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_stirling(gr_mat_t mat, int kind, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_hilbert(gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_hadamard(gr_mat_t mat, gr_ctx_t ctx);
/* todo: dft, dct */

WARN_UNUSED_RESULT int gr_mat_transpose(gr_mat_t B, const gr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_tril_classical(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_tril_recursive(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_tril_generic(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_triu_classical(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_triu_recursive(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_triu_generic(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_trace(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_mat_charpoly_berkowitz(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_charpoly_berkowitz(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_mat_charpoly_danilevsky_inplace(gr_ptr res, gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_mat_charpoly_danilevsky(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_charpoly_danilevsky(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_mat_charpoly_faddeev(gr_ptr res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_charpoly_faddeev(gr_poly_t res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_mat_charpoly_faddeev_bsgs(gr_ptr res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_charpoly_faddeev_bsgs(gr_poly_t res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_mat_charpoly_from_hessenberg(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_charpoly_from_hessenberg(gr_poly_t cp, const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_mat_charpoly_gauss(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_charpoly_gauss(gr_poly_t cp, const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_mat_charpoly_householder(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_charpoly_householder(gr_poly_t cp, const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_mat_charpoly(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_charpoly(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_hessenberg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_hessenberg_gauss(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_hessenberg_householder(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_hessenberg(const gr_mat_t mat, gr_ctx_t ctx);

int gr_mat_reduce_row(slong * column, gr_mat_t A, slong * P, slong * L, slong m, gr_ctx_t ctx);
int gr_mat_apply_row_similarity(gr_mat_t A, slong r, gr_ptr d, gr_ctx_t ctx);
int gr_mat_minpoly_field(gr_poly_t p, const gr_mat_t X, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_eigenvalues(gr_vec_t lambda, gr_vec_t mult, const gr_mat_t mat, int flags, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_eigenvalues_other(gr_vec_t lambda, gr_vec_t mult, const gr_mat_t mat, gr_ctx_t mat_ctx, int flags, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_diagonalization_precomp(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, const gr_vec_t eigenvalues, const gr_vec_t mult, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_diagonalization_generic(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, int flags, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_diagonalization(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, int flags, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_set_jordan_blocks(gr_mat_t mat, const gr_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_jordan_blocks(gr_vec_t lambda, slong * num_blocks, slong * block_lambda, slong * block_size, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_jordan_transformation(gr_mat_t mat, const gr_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_jordan_form(gr_mat_t J, gr_mat_t P, const gr_mat_t A, gr_ctx_t ctx);

truth_t gr_mat_is_scalar(const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_diagonal(const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_lower_triangular(const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_upper_triangular(const gr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_mul_diag(gr_mat_t C, const gr_mat_t A, const gr_vec_t D, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_diag_mul(gr_mat_t C, const gr_vec_t D, const gr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_exp_jordan(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_exp(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_log_jordan(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_log(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_mat_norm_max(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_norm_1(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_norm_inf(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_mat_norm_frobenius(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);

/* Test functions */

void gr_mat_test_mul(gr_method_mat_binary_op mul_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx);
void gr_mat_test_lu(gr_method_mat_lu_op lu_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx);
void gr_mat_test_det(gr_method_mat_unary_op_get_scalar det_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx);
void gr_mat_test_nonsingular_solve_tril(gr_method_mat_binary_op_with_flag solve_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx);
void gr_mat_test_nonsingular_solve_triu(gr_method_mat_binary_op_with_flag solve_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx);
void gr_mat_test_approx_mul_max_norm(gr_method_mat_binary_op mul_impl, gr_srcptr rel_tol, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx);
void gr_mat_test_approx_mul_pos_entrywise_accurate(gr_method_mat_binary_op mul_impl, gr_srcptr rel_tol, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
