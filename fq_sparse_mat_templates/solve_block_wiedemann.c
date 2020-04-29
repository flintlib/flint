/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include <string.h>
#include "templates.h"

/* Compute S_i=(M^j Y)_{0...b-1}^T for i = 0,...,ns-1 */
static void make_block_sequences(TEMPLATE(T, mat_struct)  *S, slong ns, const TEMPLATE(T, sparse_mat_t)  M, TEMPLATE(T, mat_struct)  Y[2], const TEMPLATE(T, ctx_t) ctx)
{
    slong iter, i, b = Y->c;
    TEMPLATE(T, mat_struct)  W[2];
    for (i = 0; i < 2; ++i) TEMPLATE(T, mat_window_init) (&W[i], &Y[i], 0, 0, b, b, ctx);
    for (i = iter = 0; iter < ns; ++iter, i = 1-i) 
    {
        if (iter > 0) TEMPLATE(T, sparse_mat_mul_mat) (&Y[i], M, &Y[1-i], ctx); 
        TEMPLATE(T, mat_transpose) (&S[iter], &W[i], ctx);
    }
    for (i = 0; i < 2; ++i) TEMPLATE(T, mat_window_clear) (&W[i], ctx);
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
static void coppersmith_aux_gauss(TEMPLATE(T, mat_t)  M, slong *d, const TEMPLATE(T, ctx_t) ctx) 
{
    const slong b = M->r/2;
    slong pr, pc, r, tmp;
    slong *gamma;
    TEMPLATE(T, t) cinv, cc;

    TEMPLATE(T, init) (cinv, ctx);
    TEMPLATE(T, init) (cc, ctx);

    /* Keep track of viable rows */
    gamma = flint_malloc(b*sizeof(*gamma));
    for (r = 0; r < b; ++r) gamma[r] = 1;

    for (pc = 0; pc < b; ++pc)
    {
        /* Set the pivot row to be the minimum degree row incident on column pc */
        pr = b + pc;
        for (r = 0; r < b; r++)
            if (gamma[r] && !TEMPLATE(T, is_zero) (&M->rows[r][pc], ctx) && d[r] < d[pr]) pr = r;
        if (TEMPLATE(T, is_zero) (&M->rows[pr][pc], ctx)) continue;


        /* Try to move pivot row to appropriate position (if not already there) */
        if (pr != b + pc)
        {
            tmp = d[pr]; d[pr] = d[b+pc]; d[b+pc] = tmp;

            if (!TEMPLATE(T, is_zero) (&M->rows[b + pc][pr], ctx))
                TEMPLATE(T, mat_swap_rows) (M, NULL, pr, b + pc, ctx), pr = b + pc;
            else /* Need to make new auxiliary vector and remove r from use */
                _TEMPLATE(T, vec_add) (M->rows[b + pc], M->rows[b + pc], M->rows[pr], 3*b, ctx), gamma[pr] = 0;
        }
        TEMPLATE(T, inv) (cinv, &M->rows[pr][pc], ctx);

        /* Do Gaussian elimination on first b rows */
        for (r = 0; r < b; ++r)
            if (gamma[r] && !TEMPLATE(T, is_zero) (&M->rows[r][pc], ctx))
            {
                TEMPLATE(T, mul) (cc, &M->rows[r][pc], cinv, ctx);
                TEMPLATE(T, neg) (cc, cc, ctx);
                _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T)) (M->rows[r], M->rows[pr], M->c, cc, ctx);
            }
    }
    TEMPLATE(T, clear) (cc, ctx);
    TEMPLATE(T, clear) (cinv, ctx);
    flint_free(gamma);
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
static int find_block_min_poly(TEMPLATE(T, mat_struct)  *S, slong *d, slong n, slong delta, const TEMPLATE(T, ctx_t) ctx) 
{
    int ret;
    slong t;
    slong i, k, r, b = S->r;
    slong f_len;
    TEMPLATE(T, mat_struct)  *F;
    TEMPLATE(T, mat_t)  M, D, tau, tmp;

    f_len = 1;
    F = flint_malloc((n+1)*sizeof(*F));
    TEMPLATE(T, mat_init) (&F[0], 2*b, b, ctx);
    TEMPLATE(T, mat_init) (tmp, b, b, ctx);
    for (i = 0; i < b; ++i) d[i] = 0, d[b + i] = 1, TEMPLATE(T, one) (&F[0].rows[i][i], ctx);

    /* [ D | I ] -> [ ? | tau ]*/
    TEMPLATE(T, mat_init) (M, 2*b, 3*b, ctx);

    for (t = 0, ret = -1; t < n && ret == -1; ++t)
    {
        /* Compute discrepancy matrix and tau */
        TEMPLATE(T, mat_window_init) (D, M, 0, 0, 2*b, b, ctx);
        TEMPLATE(T, mat_window_init) (tau, M, 0, b, 2*b, 3*b, ctx);
        TEMPLATE(T, mat_zero) (D, ctx);
        for (k = 0; k <= t; ++k) TEMPLATE(T, mat_addmul) (D, D, &F[k], &S[t-k], ctx);
        TEMPLATE(T, mat_one) (tau, ctx);
        TEMPLATE(T, mat_window_clear) (D, ctx);
        TEMPLATE(T, mat_window_clear) (tau, ctx);
        coppersmith_aux_gauss(M, d, ctx);

        /* Multiply F by tau * diag(I xI) */
        TEMPLATE(T, mat_window_init) (tau, M, 0, b, 2*b, 3*b, ctx); /* Needed since gauss reorders rows */
        TEMPLATE(T, mat_init) (&F[f_len++], 2*b, b, ctx);
        for (k = f_len-1; k > 0; --k)
            TEMPLATE(T, mat_mul) (&F[k], tau, &F[k-1], ctx); /* Every row multiplied by x */
        for (k = 0; k < f_len; ++k)
            for (r = 0; r < b; ++r) /* Divide first b rows by x */
            {
                if (k < f_len - 1) _TEMPLATE(T, vec_set) (F[k].rows[r], F[k+1].rows[r], b, ctx);
                else _TEMPLATE(T, vec_zero) (F[k].rows[r], b, ctx);
            }
        for (r = b; r < 2*b; ++r) _TEMPLATE(T, vec_zero) (F[0].rows[r], b, ctx), d[r] += 1;
        TEMPLATE(T, mat_window_clear) (tau, ctx);
        ret = coppersmith_stopping_criterion(d, delta, b);
    }

    /* Copy C to S, with each row reversed according to its degree */
    for (r = 0; r < b; ++r)
        for (k = 0; k <= d[r]; k++)
            _TEMPLATE(T, vec_set) (S[k].rows[r], F[d[r]-k].rows[r], b, ctx);

    for (k = 0; k < f_len; ++k) TEMPLATE(T, mat_clear) (&F[k], ctx);
    TEMPLATE(T, mat_clear) (M, ctx);
    flint_free(F);
    return ret;
}

static void make_block_sum(TEMPLATE(T, struct) *x, const TEMPLATE(T, mat_struct)  *S, const slong *d, const TEMPLATE(T, sparse_mat_t)  M, TEMPLATE(T, mat_struct)  Z[2], slong l, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, iter, b = S->r;
    slong dd;
    TEMPLATE(T, struct) *xi;

    /* Compute differences between nominal and real degree */
    dd = 0;
    while (_TEMPLATE(T, vec_is_zero) (S[dd].rows[l], b, ctx)) ++dd;

    /* Simulaneously apply all polynomials in row l to iteration of M on Z */
    xi = _TEMPLATE(T, vec_init) (M->c, ctx);
    _TEMPLATE(T, vec_zero) (x, M->c, ctx);
    for (i = iter = 0; iter <= d[l]; ++iter, i = 1 - i)
    {
        if (iter > 0) TEMPLATE(T, sparse_mat_mul_mat) (&Z[i], M, &Z[1-i], ctx);
        TEMPLATE(T, mat_mul_vec) (xi, &Z[i], S[dd + iter].rows[l], ctx);
        _TEMPLATE(T, vec_add) (x, x, xi, M->c, ctx);
    }
    _TEMPLATE(T, vec_clear) (xi, M->c, ctx);
}

int TEMPLATE(T, sparse_mat_solve_block_wiedemann) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t)  M, const TEMPLATE(T, struct) *b, slong block_size, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx)
{
    int good = 0, ret;
    TEMPLATE (T, struct) *x1;
    TEMPLATE(T, sparse_vec_t)  z;
    TEMPLATE(T, sparse_mat_t)  Mb;
    if (M->r != M->c) return 0; /* TODO */
    if (_TEMPLATE(T, vec_is_zero) (b, M->c, ctx))
    {
        _TEMPLATE(T, vec_zero) (x, M->c, ctx);
        return 1;
    }

    /* TODO: Precondition M */
    x1 = _TEMPLATE(T, vec_init) (M->c + 1, ctx);
    TEMPLATE(T, sparse_vec_init) (z, ctx);
    TEMPLATE(T, sparse_mat_init) (Mb, M->r, M->c, ctx);
    TEMPLATE(T, sparse_mat_set) (Mb, M, ctx);
    TEMPLATE(T, sparse_mat_append_col) (Mb, b, ctx);
    TEMPLATE(T, sparse_mat_append_row) (Mb, z, ctx);

    ret = TEMPLATE(T, sparse_mat_nullvector_block_wiedemann) (x1, Mb, block_size, state, ctx);
    if (ret && !TEMPLATE(T, is_zero) (&x1[M->c], ctx)) 
    {
        TEMPLATE(T, inv) (&x1[M->c], &x1[M->c], ctx);
        TEMPLATE(T, neg) (&x1[M->c], &x1[M->c], ctx);
        _TEMPLATE(T, TEMPLATE(vec_scalar_mul, T)) (x, x1, M->c, &x1[M->c], ctx);
        good = 1;
    }
    TEMPLATE(T, sparse_vec_clear) (z, ctx);
    TEMPLATE(T, sparse_mat_clear) (Mb, ctx);
    _TEMPLATE(T, vec_clear) (x1, M->c + 1, ctx);
    return good;
}

int TEMPLATE(T, sparse_mat_nullvector_block_wiedemann) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t)  M, slong block_size, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx) 
{
    int ret = 0;
    slong l, ns, k;
    slong *d;
    TEMPLATE(T, struct) *b;
    TEMPLATE(T, mat_struct)  Y[3], *S;
    if (M->r != M->c) return 0; /* TODO */

    ns = 2*M->r/block_size + 3; /* Maybe 5? */
    S = flint_malloc(ns*sizeof(*S));
    d = flint_calloc(2*block_size, sizeof(*d));
    b = _TEMPLATE(T, vec_init) (M->r, ctx);
    for (k = 0; k < ns; ++k) TEMPLATE(T, mat_init) (&S[k], block_size, block_size, ctx);
    for (l = 0; l < 3; ++l) TEMPLATE(T, mat_init) (&Y[l], M->c, block_size, ctx);
    do TEMPLATE(T, mat_randtest) (&Y[0], state, ctx);
    while (TEMPLATE(T, mat_is_zero) (&Y[0], ctx));

    TEMPLATE(T, sparse_mat_mul_mat) (&Y[1], M, &Y[0], ctx);
    make_block_sequences(S, ns, M, &Y[1], ctx);
    find_block_min_poly(S, d, ns, M->r, ctx);

    for (l = 0; l < block_size; ++l)
    {
        TEMPLATE(T, mat_set) (&Y[1], &Y[0], ctx);
        make_block_sum(x, S, d, M, Y + 1, l, ctx);
        TEMPLATE(T, sparse_mat_mul_vec) (b, M, x, ctx);
        if (!_TEMPLATE(T, vec_is_zero) (x, M->c, ctx) && _TEMPLATE(T, vec_is_zero) (b, M->r, ctx)) {ret = 1; break;};
    }    
    return ret;
}

#endif
