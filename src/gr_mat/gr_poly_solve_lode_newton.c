/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"
#include "gr_poly.h"

/*
 * This is a FLINT implementation of the algorithm described in:
 * Fast computation of power series solutions of systems of differential equations
 * by A. Bostan, F. Chyzak, F. Ollivier, B. Salvy, É. Schost and A. Sedoglavic
 * in Proceedings of the Eighteenth Annual ACM-SIAM Symposium on Discrete Algorithms, 2007, pp. 1012 -- 1021.
 *
 * Compared with the algorithm cited above, we use an iteration for the inverse of the denominator of the matrix A.
 * When A is a companion matrix, the polynomial matrix multiplication by A modulo z^m is optimized.
 * An optimal sequence of precisions is used to reach the target length, even when it is not a power of two.
 */

/*
 * TODO: Newton iteration context struct? For keeping Z and A_denominator_inv.
 */

/* This is a specialized implementation of univariate polynomial matrix multiplication
 * for the case where A = (companion matrix) and the result is computed modulo z^m.
 * [ 0 1 0 0 ]   [ b00 b01 b02 b03 ]   [ b10 b11 b12 b13 ]
 * [ 0 0 1 0 ] * [ b10 b11 b12 b13 ] = [ b20 b21 b22 b23 ]
 * [ 0 0 0 1 ]   [ b20 b21 b22 b23 ]   [ b30 b31 b32 b33 ]
 * [ a b c d ]   [ b30 b31 b32 b33 ]   [ ... ... ... ... ]
 */
static int
gr_mat_gr_poly_mullow_companion(gr_mat_t res, gr_mat_t A, gr_mat_t B, slong len, gr_ctx_t poly_ctx)
{
    gr_ctx_struct *coeff_ctx = POLYNOMIAL_CTX(poly_ctx)->base_ring;
    gr_mat_t tmp_mat;
    gr_poly_t tmp_poly;
    int status = GR_SUCCESS;

    gr_mat_init(tmp_mat, A->r, A->c, poly_ctx);
    gr_poly_init(tmp_poly, coeff_ctx);

    status |= gr_mat_zero(tmp_mat, poly_ctx);
    for (int i = 0; i < A->r - 1; i++)
    {
        for (int j = 0; j < A->c; j++)
        {
            gr_srcptr entry_B = gr_mat_entry_srcptr(B, i + 1, j, poly_ctx);
            gr_ptr entry_tmp = gr_mat_entry_ptr(tmp_mat, i, j, poly_ctx);
            status |= gr_poly_set(entry_tmp, entry_B, coeff_ctx);
        }
    }

    for (int j = 0; j < A->c; j++)
    {
        gr_ptr entry_tmp = gr_mat_entry_ptr(tmp_mat, A->r - 1, j, poly_ctx);
        for (int i = 0; i < B->r; i++)
        {
            gr_srcptr entry_A = gr_mat_entry_srcptr(A, A->r - 1, i, poly_ctx);
            gr_srcptr entry_B = gr_mat_entry_srcptr(B, i, j, poly_ctx);
            status |= gr_poly_mullow(tmp_poly, entry_A, entry_B, len, coeff_ctx);
            status |= gr_poly_add(entry_tmp, entry_tmp, tmp_poly, coeff_ctx);
        }
    }

    gr_mat_swap(res, tmp_mat, poly_ctx);

    gr_poly_clear(tmp_poly, coeff_ctx);
    gr_mat_clear(tmp_mat, poly_ctx);
    return status;
}

int
_gr_mat_gr_poly_solve_lode_newton_start(gr_mat_t Y, gr_mat_t Z, gr_poly_t A_denominator_inv, const gr_mat_t A_numerator, const gr_poly_t A_denominator, const gr_mat_t Y0, gr_ctx_t sol_poly_ctx)
{
    gr_ctx_struct *sol_coeff_ctx = POLYNOMIAL_CTX(sol_poly_ctx)->base_ring;

    gr_poly_t tmp_poly;
    gr_mat_t Y0inv, tmp_mat;
    slong i, j;
    int m, status = GR_SUCCESS;

    /* Calculate the first term of the power series inverse of A_denominator: 1 / A_denominator(0). */
    status |= gr_poly_one(A_denominator_inv, sol_coeff_ctx);
    if (gr_poly_length(A_denominator, sol_coeff_ctx) > 0)
        status |= gr_poly_div_scalar(A_denominator_inv, A_denominator_inv, gr_poly_entry_srcptr(A_denominator, 0, sol_coeff_ctx), sol_coeff_ctx);
    else
        status = GR_DOMAIN;

    if (status != GR_SUCCESS)
        return status;

    /* Check that Y0 is invertible. */
    gr_mat_init(Y0inv, Y0->r, Y0->c, sol_coeff_ctx);
    status |= gr_mat_inv(Y0inv, Y0, sol_coeff_ctx);
    if (status != GR_SUCCESS)
    {
        gr_mat_clear(Y0inv, sol_coeff_ctx);
        return status;
    }

    gr_mat_init(tmp_mat, Y->r, Y->c, sol_poly_ctx);
    gr_poly_init(tmp_poly, sol_coeff_ctx);

    /* Z = Y0^{-1}; this will become the inverse of Y as a power series. */
    status |= gr_mat_set_gr_mat_other(Z, Y0inv, sol_coeff_ctx, sol_poly_ctx);

    /* Next we are going to set Y = (I + A(0)*z)*Y0 in a few steps. */

    /* Y = A(0)*z */
    for (i = 0; i < Y->r; i++)
    {
        for (j = 0; j < Y->c; j++)
        {
            gr_srcptr entry_A = gr_mat_entry_srcptr(A_numerator, i, j, sol_poly_ctx);
            gr_ptr entry_Y = gr_mat_entry_ptr(Y, i, j, sol_poly_ctx);
            status |= gr_poly_truncate(entry_Y, entry_A, 1, sol_coeff_ctx);
            status |= gr_poly_shift_left(entry_Y, entry_Y, 1, sol_coeff_ctx);
        }
    }
    status |= gr_mat_mul_scalar(Y, Y, A_denominator_inv, sol_poly_ctx);

    /* Y = (I + Y)*Y0 */
    status |= gr_mat_add_ui(Y, Y, 1, sol_poly_ctx);
    status |= gr_mat_set_gr_mat_other(tmp_mat, Y0, sol_coeff_ctx, sol_poly_ctx);
    status |= gr_mat_mul(Y, Y, tmp_mat, sol_poly_ctx);

    /* One initial iteration for A_denominator_inv mod z^2. */
    m = 2;
    status |= gr_poly_truncate(tmp_poly, A_denominator, m, sol_coeff_ctx);
    status |= gr_poly_mullow(tmp_poly, tmp_poly, A_denominator_inv, m, sol_coeff_ctx);
    status |= gr_poly_shift_right(tmp_poly, tmp_poly, m / 2, sol_coeff_ctx);
    status |= gr_poly_mullow(tmp_poly, tmp_poly, A_denominator_inv, m / 2, sol_coeff_ctx);
    status |= gr_poly_shift_left(tmp_poly, tmp_poly, m / 2, sol_coeff_ctx);
    status |= gr_poly_sub(A_denominator_inv, A_denominator_inv, tmp_poly, sol_coeff_ctx);

    gr_poly_clear(tmp_poly, sol_coeff_ctx);
    gr_mat_clear(tmp_mat, sol_poly_ctx);
    gr_mat_clear(Y0inv, sol_coeff_ctx);

    return status;
}

int
_gr_mat_gr_poly_solve_lode_newton_step(gr_mat_t Y, gr_mat_t Z, gr_poly_t A_denominator_inv, slong len, const gr_mat_t A_numerator, const gr_poly_t A_denominator, int A_is_companion, gr_ctx_t sol_poly_ctx)
{
    gr_ctx_struct *sol_coeff_ctx = POLYNOMIAL_CTX(sol_poly_ctx)->base_ring;

    gr_mat_t Err, tmp_mat;
    gr_poly_t tmp_poly;
    slong i, j, m = (len + 1) / 2;
    int status = GR_SUCCESS;

    gr_mat_init(Err, Y->r, Y->c, sol_poly_ctx);
    gr_mat_init(tmp_mat, Y->r, Y->c, sol_poly_ctx);
    gr_poly_init(tmp_poly, sol_coeff_ctx);

    /* Z = Z + (Z*(I - Y*Z) mod z^m) */
    status |= gr_mat_mul(tmp_mat, Y, Z, sol_poly_ctx);
    status |= gr_mat_sub_ui(tmp_mat, tmp_mat, 1, sol_poly_ctx);
    status |= gr_mat_mul(tmp_mat, Z, tmp_mat, sol_poly_ctx);
    for (i = 0; i < tmp_mat->r; i++)
    {
        for (j = 0; j < tmp_mat->c; j++)
        {
            gr_ptr entry_tmp = gr_mat_entry_ptr(tmp_mat, i, j, sol_poly_ctx);
            status |= gr_poly_truncate(entry_tmp, entry_tmp, m, sol_coeff_ctx);
        }
    }
    status |= gr_mat_sub(Z, Z, tmp_mat, sol_poly_ctx);

    /* Iteration for A_denominator_inv mod z^len. */
    status |= gr_poly_truncate(tmp_poly, A_denominator, len, sol_coeff_ctx);
    /* TODO: Use middle product here? */
    status |= gr_poly_mullow(tmp_poly, tmp_poly, A_denominator_inv, len, sol_coeff_ctx);
    status |= gr_poly_shift_right(tmp_poly, tmp_poly, m, sol_coeff_ctx);
    status |= gr_poly_mullow(tmp_poly, tmp_poly, A_denominator_inv, m, sol_coeff_ctx);
    status |= gr_poly_shift_left(tmp_poly, tmp_poly, m, sol_coeff_ctx);
    status |= gr_poly_sub(A_denominator_inv, A_denominator_inv, tmp_poly, sol_coeff_ctx);

    /* Next we are going to set Err = Y' - (A mod z^{len - 1})*Y mod z^{len - 1} in a few steps. */

    /* tmp_mat = (A mod z^{len - 1})*Y */
    for (i = 0; i < Y->r; i++)
    {
        for (j = 0; j < Y->c; j++)
        {
            gr_ptr entry_tmp = gr_mat_entry_ptr(tmp_mat, i, j, sol_poly_ctx);
            status |= gr_poly_truncate(entry_tmp, gr_mat_entry_srcptr(A_numerator, i, j, sol_poly_ctx), len - 1, sol_coeff_ctx);
            status |= gr_poly_mullow(entry_tmp, entry_tmp, A_denominator_inv, len - 1, sol_coeff_ctx);
        }
    }
    if (A_is_companion)
        status |= gr_mat_gr_poly_mullow_companion(tmp_mat, tmp_mat, Y, len - 1, sol_poly_ctx);
    else
        status |= gr_mat_mul(tmp_mat, tmp_mat, Y, sol_poly_ctx);

    /* Err = Y' - tmp_mat */
    for (i = 0; i < Y->r; i++)
    {
        for (j = 0; j < Y->c; j++)
        {
            gr_srcptr entry_Y = gr_mat_entry_srcptr(Y, i, j, sol_poly_ctx);
            gr_ptr entry_Err = gr_mat_entry_ptr(Err, i, j, sol_poly_ctx);
            status |= gr_poly_derivative(entry_Err, entry_Y, sol_coeff_ctx);
        }
    }
    status |= gr_mat_sub(Err, Err, tmp_mat, sol_poly_ctx);

    /* Y = Y - (Y*integral(Z*Err) mod z^len) */
    status |= gr_mat_mul(tmp_mat, Z, Err, sol_poly_ctx);
    for (i = 0; i < Y->r; i++)
    {
        for (j = 0; j < Y->c; j++)
        {
            gr_ptr entry_tmp = gr_mat_entry_ptr(tmp_mat, i, j, sol_poly_ctx);
            status |= gr_poly_integral(entry_tmp, entry_tmp, sol_coeff_ctx);
        }
    }
    status |= gr_mat_mul(tmp_mat, Y, tmp_mat, sol_poly_ctx);
    for (i = 0; i < tmp_mat->r; i++)
    {
        for (j = 0; j < tmp_mat->c; j++)
        {
            gr_ptr entry_tmp = gr_mat_entry_ptr(tmp_mat, i, j, sol_poly_ctx);
            status |= gr_poly_truncate(entry_tmp, entry_tmp, len, sol_coeff_ctx);
        }
    }
    status |= gr_mat_sub(Y, Y, tmp_mat, sol_poly_ctx);

    gr_poly_clear(tmp_poly, sol_coeff_ctx);
    gr_mat_clear(tmp_mat, sol_poly_ctx);
    gr_mat_clear(Err, sol_poly_ctx);

    return status;
}

int
_gr_mat_gr_poly_solve_lode_newton(gr_mat_t Y, gr_mat_t Z, const gr_mat_t A_numerator, const gr_poly_t A_denominator, const gr_mat_t Y0, slong len, gr_ctx_t A_poly_ctx, gr_ctx_t sol_poly_ctx)
{
    gr_ctx_struct *A_coeff_ctx = POLYNOMIAL_CTX(A_poly_ctx)->base_ring;
    gr_ctx_struct *sol_coeff_ctx = POLYNOMIAL_CTX(sol_poly_ctx)->base_ring;

    gr_mat_t A_numerator_sol;
    gr_poly_t A_denominator_sol, A_denominator_inv;
    slong i, j, n;
    slong a[FLINT_BITS];
    int A_is_companion = 1, status = GR_SUCCESS;

    /* Detect companion case. */
    for (i = 0; i < A_numerator->r - 1; i++)
    {
        for (j = 0; j < A_numerator->c; j++)
        {
            gr_srcptr entry_A = gr_mat_entry_srcptr(A_numerator, i, j, A_poly_ctx);
            if (j == i + 1)
                A_is_companion &= (gr_poly_equal(entry_A, A_denominator, A_coeff_ctx) == T_TRUE);
            else
                A_is_companion &= (gr_poly_is_zero(entry_A, A_coeff_ctx) == T_TRUE);
        }
    }

    gr_poly_init(A_denominator_sol, sol_coeff_ctx);
    gr_poly_init(A_denominator_inv, sol_coeff_ctx);

    /* Convert A_numerator and A_denominator into the solution context. */
    if (A_poly_ctx == sol_poly_ctx)
    {
        /* TODO: Avoid copying? */
        status |= gr_poly_set(A_denominator_sol, A_denominator, sol_coeff_ctx);
        status |= gr_mat_init_set(A_numerator_sol, A_numerator, sol_poly_ctx);
    }
    else
    {
        status |= gr_poly_set_gr_poly_other(A_denominator_sol, A_denominator, A_coeff_ctx, sol_coeff_ctx);
        gr_mat_init(A_numerator_sol, A_numerator->r, A_numerator->c, sol_poly_ctx);
        status |= gr_mat_set_gr_mat_other(A_numerator_sol, A_numerator, A_poly_ctx, sol_poly_ctx);
    }

    /* TODO: Add basecase? */

    a[i = 0] = n = len;
    while (n > 1)
        a[++i] = (n = (n + 1) / 2);

    status |= _gr_mat_gr_poly_solve_lode_newton_start(Y, Z, A_denominator_inv, A_numerator_sol, A_denominator_sol, Y0, sol_poly_ctx);
    if (status == GR_SUCCESS)
        for (i--; i >= 0; i--)
            status |= _gr_mat_gr_poly_solve_lode_newton_step(Y, Z, A_denominator_inv, a[i], A_numerator_sol, A_denominator_sol, A_is_companion, sol_poly_ctx);

    /* TODO: A way to save A_denominator_inv? */

    gr_mat_clear(A_numerator_sol, sol_poly_ctx);
    gr_poly_clear(A_denominator_inv, sol_coeff_ctx);
    gr_poly_clear(A_denominator_sol, sol_coeff_ctx);

    return status;
}

int
gr_mat_gr_poly_solve_lode_newton(gr_mat_t Y, const gr_mat_t A_numerator, const gr_poly_t A_denominator, const gr_mat_t Y0, slong len, gr_ctx_t A_poly_ctx, gr_ctx_t sol_poly_ctx)
{
    gr_mat_t Z;
    int status;

    gr_mat_init(Z, Y0->r, Y0->c, sol_poly_ctx);
    status = _gr_mat_gr_poly_solve_lode_newton(Y, Z, A_numerator, A_denominator, Y0, len, A_poly_ctx, sol_poly_ctx);
    gr_mat_clear(Z, sol_poly_ctx);

    return status;
}
