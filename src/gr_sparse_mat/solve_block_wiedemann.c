/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "gr_sparse_mat.h"

/* Compute S_i=(M^j Y)_{0...b-1}^T for i = 0,...,ns-1 */
static int make_block_sequences(gr_mat_struct *S, slong ns, const gr_lil_mat_t M, gr_mat_struct Y[2], gr_ctx_t ctx) 
{
    slong iter, i, b = Y->c;
    gr_mat_struct W[2];
    int status = GR_SUCCESS;

    for (i = 0; i < 2; ++i) gr_mat_window_init(&W[i], &Y[i], 0, 0, b, b, ctx);
    for (i = iter = 0; iter < ns; ++iter, i = 1-i) 
    {
        if (iter > 0) status |= gr_lil_mat_mul_mat(&Y[i], M, &Y[1-i], ctx); 
        status |= gr_mat_transpose(&S[iter], &W[i], ctx);
    }
    for (i = 0; i < 2; ++i) gr_mat_window_clear(&W[i], ctx);
    return status;
}

/**
 * Run Guassian elimination on the first b columns of the augmented 
 * matrix M = [ D | I], yielding a final matrix 
 *              [   |       ]    [ Z |       ]
 *              [ D |   I   ] -> [---|  tau  ]
 *              [   |       ]    [ L |       ]
 * where the the number of nonzero rows in Z is the ith rank. We choose 
 * the pivot row for a given column to be the one with minimal degree.
**/
static int coppersmith_aux_gauss(gr_mat_t M, slong *d, gr_ctx_t ctx) 
{
    const slong b = M->r/2;
    slong pr, pc, r, tmp;
    slong *gamma;
    gr_ptr cinv, coeff;
    int status = GR_SUCCESS;

    /* Keep track of viable rows */
    gamma = flint_malloc(b * sizeof(slong));
    for (r = 0; r < b; ++r) gamma[r] = 1;

    GR_TMP_INIT2(cinv, coeff, ctx);
    for (pc = 0; pc < b; ++pc)
    {
        /* Set the pivot row to be the minimum degree row incident on column pc */
        pr = b + pc;
        for (r = 0; r < b; r++)
            if (gamma[r] && gr_mat_entry_ptr(M, r, pc, ctx) && d[r] < d[pr]) pr = r;
        if (gr_is_zero(gr_mat_entry_ptr(M, pr, pc, ctx), ctx) == T_TRUE) continue;


        /* Try to move pivot row to appropriate position (if not already there) */
        if (pr != b + pc)
        {
            tmp = d[pr]; d[pr] = d[b+pc]; d[b+pc] = tmp;

            if (gr_mat_entry_ptr(M, b + pc, pr, ctx)) 
                status |= gr_mat_swap_rows(M, NULL, pr, b + pc, ctx), pr = b + pc;
            else /* Need to make new auxiliary vector and remove r from use */
                status |= _gr_vec_add(M->rows[b + pc], M->rows[b + pc], M->rows[pr], 3*b, ctx), gamma[pr] = 0;
        }
        status |= gr_inv(cinv, gr_mat_entry_ptr(M, pr, pc, ctx), ctx);

        /* Do Gaussian elimination on first b rows */
        for (r = 0; r < b; ++r)
            if (gamma[r] && gr_is_zero(gr_mat_entry_ptr(M, r, pc, ctx), ctx) == T_FALSE)
            {
                status |= gr_mul(coeff, gr_mat_entry_ptr(M, r, pc, ctx), cinv, ctx);
                status |= gr_neg(coeff, coeff, ctx);
                status |= _gr_vec_addmul_scalar(
                    M->rows[r], M->rows[pr], M->c,
                    coeff, ctx
                );

            }
    }
    flint_free(gamma);
    GR_TMP_CLEAR2(cinv, coeff, ctx);
    return status;
}

/* Stop with failure if sum(d_0 ... d_{b-1}) < delta */
/* Stop with success if sum(d_0 ... d_{b-1}) < delta + max(d_0 ... d_{b-1}) - min(d_b ... d_{2b-1})  */
static int coppersmith_stopping_criterion(slong *d, slong delta, slong b)
{
    slong tmp, r;

    /* Sum degrees of generating polynomials */
    tmp = d[0]; for (r = 1; r < b; ++r) tmp += d[r];
    delta -= tmp;
    if (delta < 0) return 0; /* Insufficient degree */

    /* Add maximum degree of first b polys and subtract minimum degree of last b */
    tmp = d[0]; for (r = 1; r < b; ++r) if (d[r] > tmp) tmp = d[r];
    delta += tmp;
    tmp = d[b]; for (r = b + 1; r < 2*b; ++r) if (d[r] < tmp) tmp = d[r];
    delta -= tmp;
    return delta < 0 ? 1 : -1;
}

/**
 * Generalization of Berlekamp-Massey due to Coppersmith.
 * Iteratively computes a sequence F representing 2b polynomials:
 *  - the first b are the current (reversed) generating polynomials
 *  - the last b are certain auxiliary polynomials.
**/
static int find_block_min_poly(gr_mat_struct *S, slong *d, slong n, slong delta, gr_ctx_t ctx) 
{
    int ret;
    slong t;
    slong i, k, r, b = S->r;
    slong f_len;
    gr_mat_struct *F;
    gr_mat_t M, D, tau, tmp;
    int status = GR_SUCCESS;

    f_len = 1;
    F = flint_malloc((n + 1) * sizeof(gr_mat_struct));
    gr_mat_init(&F[0], 2*b, b, ctx);
    gr_mat_init(tmp, b, b, ctx);
    for (i = 0; i < b; ++i)
    {
        d[i] = 0;
        d[b + i] = 1;
        status |= gr_one(gr_mat_entry_ptr(&F[0], i, i, ctx), ctx);
    }

    /* [ D | I ] -> [ ? | tau ]*/
    gr_mat_init(M, 2*b, 3*b, ctx);

    for (t = 0, ret = -1; t < n && ret == -1; ++t)
    {
        /* Compute discrepancy matrix and tau */
        gr_mat_window_init(D, M, 0, 0, 2*b, b, ctx);
        gr_mat_window_init(tau, M, 0, b, 2*b, 3*b, ctx);
        status |= gr_mat_zero(D, ctx);
        for (k = 0; k <= t; ++k) status |= gr_mat_addmul(D, D, &F[k], &S[t-k], ctx);
        status |= gr_mat_one(tau, ctx);
        gr_mat_window_clear(D, ctx);
        gr_mat_window_clear(tau, ctx);
        status |= coppersmith_aux_gauss(M, d, ctx);

        /* Multiply F by tau * diag(I xI) */
        gr_mat_window_init(tau, M, 0, b, 2*b, 3*b, ctx); /* Needed since gauss reorders rows */
        gr_mat_init(&F[f_len++], 2*b, b, ctx);
        for (k = f_len-1; k > 0; --k)
            status |= gr_mat_mul(&F[k], tau, &F[k-1], ctx); /* Every row multiplied by x */
        for (k = 0; k < f_len; ++k)
            for (r = 0; r < b; ++r) /* Divide first b rows by x */
            {
                if (k < f_len - 1)
                    status |= _gr_vec_set(F[k].rows[r], F[k+1].rows[r], b, ctx);
                else
                    status |= _gr_vec_zero(F[k].rows[r], b, ctx);
            }
        for (r = b; r < 2*b; ++r)
            status |= _gr_vec_zero(F[0].rows[r], b, ctx), d[r] += 1;
        gr_mat_window_clear(tau, ctx);
        ret = coppersmith_stopping_criterion(d, delta, b);
    }

    /* Copy C to S, with each row reversed according to its degree */
    for (r = 0; r < b; ++r)
        for (k = 0; k <= d[r]; k++)
            status |= _gr_vec_set(S[k].rows[r], F[d[r]-k].rows[r], b, ctx);

    for (k = 0; k < f_len; ++k)
        gr_mat_clear(&F[k], ctx);
    gr_mat_clear(M, ctx);
    flint_free(F);
    return status;
}

static int make_block_sum(gr_ptr x, const gr_mat_struct *S, const slong *d, const gr_lil_mat_t M, gr_mat_struct Z[2], slong l, gr_ctx_t ctx)
{
    slong i, iter, b = S->r;
    slong dd;
    gr_ptr xi;
    int status = GR_SUCCESS;

    /* Compute differences between nominal and real degree */
    dd = 0;
    while (_gr_vec_is_zero(S[dd].rows[l], b, ctx) == T_TRUE) ++dd;

    /* Simulaneously apply all polynomials in row l to iteration of M on Z */
    GR_TMP_INIT_VEC(xi, M->c, ctx);
    status |= _gr_vec_zero(x, M->c, ctx);
    for (i = iter = 0; iter <= d[l]; ++iter, i = 1 - i)
    {
        if (iter > 0) status |= gr_lil_mat_mul_mat(&Z[i], M, &Z[1-i], ctx);
        status |= gr_mat_mul_vec(xi, &Z[i], S[dd + iter].rows[l], ctx);
        status |= _gr_vec_add(x, x, xi, M->c, ctx);
    }
    GR_TMP_CLEAR_VEC(xi, M->c, ctx);
    return status;
}

int gr_lil_mat_solve_block_wiedemann(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, slong block_size, flint_rand_t state, gr_ctx_t ctx)
{
    int i;
    slong sz = ctx->sizeof_elem;
    gr_ptr x1, coeff;
    gr_lil_mat_t Mb;
    int status = GR_SUCCESS;

    if (M->r != M->c) return 0; /* TODO */
    if (_gr_vec_is_zero(b, M->c, ctx) == T_TRUE)
    {
        return _gr_vec_zero(x, M->c, ctx);
    }

    /* TODO: Precondition M */
    GR_TMP_INIT(coeff, ctx);
    GR_TMP_INIT_VEC(x1, M->c + 1, ctx);
    gr_lil_mat_init(Mb, M->r + 1, M->c + 1, ctx);
    for (i = 0; i < M->r; ++i)
    {
        gr_sparse_vec_fit_nnz(&Mb->rows[i], M->rows[i].nnz + 1, ctx);
        status |= gr_sparse_vec_set(&Mb->rows[i], &M->rows[i], ctx);
        status |= gr_sparse_vec_set_entry(&Mb->rows[i], M->c, GR_ENTRY(b, i, sz), ctx);
    }
    status |= gr_lil_mat_set(Mb, M, ctx);

    status |= gr_lil_mat_nullvector_block_wiedemann(x1, Mb, block_size, state, ctx);
    if (status == GR_SUCCESS && gr_is_zero(GR_ENTRY(x1, M->c, sz), ctx) == T_FALSE) 
    {
        status |= gr_inv(coeff, GR_ENTRY(x1, M->c, sz), ctx);
        status |= gr_neg(coeff, coeff, ctx);
        status |= _gr_vec_mul_scalar(x, x1, M->c, coeff, ctx);
    }
    gr_lil_mat_clear(Mb, ctx);
    GR_TMP_CLEAR_VEC(x1, M->c + 1, ctx);
    GR_TMP_CLEAR(coeff, ctx);
    return status;
}

int gr_lil_mat_nullvector_block_wiedemann(gr_ptr x, const gr_lil_mat_t M, slong block_size, flint_rand_t state, gr_ctx_t ctx) 
{
    slong l, ns, k;
    slong *d;
    gr_ptr b;
    gr_mat_struct Y[3], *S;
    int status = GR_SUCCESS;
    if (M->r != M->c) return 0; /* TODO */

    ns = 2*M->r/block_size + 3; /* Maybe 5? */
    S = flint_malloc(ns*sizeof(*S));
    d = flint_calloc(2*block_size, sizeof(*d));
    GR_TMP_INIT_VEC(b, M->r, ctx);
    for (k = 0; k < ns; ++k)
        gr_mat_init(&S[k], block_size, block_size, ctx);
    for (l = 0; l < 3; ++l)
        gr_mat_init(&Y[l], M->c, block_size, ctx);
    do
        status |= gr_mat_randtest(&Y[0], state, ctx);
    while (gr_mat_is_zero(&Y[0], ctx) == T_TRUE);

    status |= gr_lil_mat_mul_mat(&Y[1], M, &Y[0], ctx);
    status |= make_block_sequences(S, ns, M, &Y[1], ctx);
    status |= find_block_min_poly(S, d, ns, M->r, ctx);

    for (l = 0; l < block_size; ++l)
    {
        status |= gr_mat_set(&Y[1], &Y[0], ctx);
        status |= make_block_sum(x, S, d, M, Y + 1, l, ctx);
        status |= gr_lil_mat_mul_vec(b, M, x, ctx);
        if (!_gr_vec_is_zero(x, M->c, ctx) && _gr_vec_is_zero(b, M->r, ctx)) {status = GR_DOMAIN; break;};
    }    
    GR_TMP_CLEAR_VEC(b, M->r, ctx);
    return status;
}


