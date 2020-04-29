/*
    Copyright (C) 2010 Fredrik Johansson

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
#include "nmod_sparse_mat.h"

/* Compute S_i=(M^j Y)_{0...b-1}^T for i = 0,...,ns-1 */
static void make_block_sequences(nmod_mat_struct *S, slong ns, const nmod_sparse_mat_t M, nmod_mat_struct Y[2]) 
{
    slong iter, i, b = Y->c;
    nmod_mat_struct W[2];
    for (i = 0; i < 2; ++i) nmod_mat_window_init(&W[i], &Y[i], 0, 0, b, b);
    for (i = iter = 0; iter < ns; ++iter, i = 1-i) 
    {
        if (iter > 0) nmod_sparse_mat_mul_mat(&Y[i], M, &Y[1-i]); 
        nmod_mat_transpose(&S[iter], &W[i]);
    }
    for (i = 0; i < 2; ++i) nmod_mat_window_clear(&W[i]);
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
static void coppersmith_aux_gauss(nmod_mat_t M, slong *d) 
{
    const slong b = M->r/2;
    slong pr, pc, r, tmp;
    slong *gamma;
    mp_limb_t cinv;

    /* Keep track of viable rows */
    gamma = flint_malloc(b*sizeof(*gamma));
    for (r = 0; r < b; ++r) gamma[r] = 1;

    for (pc = 0; pc < b; ++pc)
    {
        /* Set the pivot row to be the minimum degree row incident on column pc */
        pr = b + pc;
        for (r = 0; r < b; r++)
            if (gamma[r] && M->rows[r][pc] && d[r] < d[pr]) pr = r;
        if (M->rows[pr][pc] == UWORD(0)) continue;


        /* Try to move pivot row to appropriate position (if not already there) */
        if (pr != b + pc)
        {
            tmp = d[pr]; d[pr] = d[b+pc]; d[b+pc] = tmp;

            if (M->rows[b + pc][pr]) 
                nmod_mat_swap_rows(M, NULL, pr, b + pc), pr = b + pc;
            else /* Need to make new auxiliary vector and remove r from use */
                _nmod_vec_add(M->rows[b + pc], M->rows[b + pc], M->rows[pr], 3*b, M->mod), gamma[pr] = 0;
        }
        cinv = nmod_inv(M->rows[pr][pc], M->mod);

        /* Do Gaussian elimination on first b rows */
        for (r = 0; r < b; ++r)
            if (gamma[r] && M->rows[r][pc])
                _nmod_vec_scalar_addmul_nmod(M->rows[r], M->rows[pr], M->c,
                                             nmod_neg(nmod_mul(M->rows[r][pc], cinv, M->mod), M->mod), M->mod);
    }
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
static int find_block_min_poly(nmod_mat_struct *S, slong *d, slong n, slong delta) 
{
    int ret;
    slong t;
    slong i, k, r, b = S->r;
    slong f_len;
    nmod_mat_struct *F;
    nmod_mat_t M, D, tau, tmp;

    f_len = 1;
    F = flint_malloc((n+1)*sizeof(*F));
    nmod_mat_init(&F[0], 2*b, b, S->mod.n);
    nmod_mat_init(tmp, b, b, S->mod.n);
    for (i = 0; i < b; ++i) d[i] = 0, d[b + i] = 1, F[0].rows[i][i] = 1;

    /* [ D | I ] -> [ ? | tau ]*/
    nmod_mat_init(M, 2*b, 3*b, S->mod.n);

    for (t = 0, ret = -1; t < n && ret == -1; ++t)
    {
        /* Compute discrepancy matrix and tau */
        nmod_mat_window_init(D, M, 0, 0, 2*b, b);
        nmod_mat_window_init(tau, M, 0, b, 2*b, 3*b);
        nmod_mat_zero(D);
        for (k = 0; k <= t; ++k) nmod_mat_addmul(D, D, &F[k], &S[t-k]);
        nmod_mat_one(tau);
        nmod_mat_window_clear(D);
        nmod_mat_window_clear(tau);
        coppersmith_aux_gauss(M, d);

        /* Multiply F by tau * diag(I xI) */
        nmod_mat_window_init(tau, M, 0, b, 2*b, 3*b); /* Needed since gauss reorders rows */
        nmod_mat_init(&F[f_len++], 2*b, b, S->mod.n);
        for (k = f_len-1; k > 0; --k)
            nmod_mat_mul(&F[k], tau, &F[k-1]); /* Every row multiplied by x */
        for (k = 0; k < f_len; ++k)
            for (r = 0; r < b; ++r) /* Divide first b rows by x */
            {
                if (k < f_len - 1) _nmod_vec_set(F[k].rows[r], F[k+1].rows[r], b);
                else _nmod_vec_zero(F[k].rows[r], b);
            }
        for (r = b; r < 2*b; ++r) _nmod_vec_zero(F[0].rows[r], b), d[r] += 1;
        nmod_mat_window_clear(tau);
        ret = coppersmith_stopping_criterion(d, delta, b);
    }

    /* Copy C to S, with each row reversed according to its degree */
    for (r = 0; r < b; ++r)
        for (k = 0; k <= d[r]; k++)
            _nmod_vec_set(S[k].rows[r], F[d[r]-k].rows[r], b);

    for (k = 0; k < f_len; ++k) nmod_mat_clear(&F[k]);
    nmod_mat_clear(M);
    flint_free(F);
    return ret;
}

static void make_block_sum(mp_ptr x, const nmod_mat_struct *S, const slong *d, const nmod_sparse_mat_t M, nmod_mat_struct Z[2], slong l)
{
    slong i, iter, b = S->r;
    slong dd;
    mp_ptr xi;

    /* Compute differences between nominal and real degree */
    dd = 0;
    while (_nmod_vec_is_zero(S[dd].rows[l], b)) ++dd;

    /* Simulaneously apply all polynomials in row l to iteration of M on Z */
    xi = _nmod_vec_init(M->c);
    _nmod_vec_zero(x, M->c);
    for (i = iter = 0; iter <= d[l]; ++iter, i = 1 - i)
    {
        if (iter > 0) nmod_sparse_mat_mul_mat(&Z[i], M, &Z[1-i]);
        nmod_mat_mul_vec(xi, &Z[i], S[dd + iter].rows[l]);
        _nmod_vec_add(x, x, xi, M->c, M->mod);
    }
    _nmod_vec_clear(xi);
}

int nmod_sparse_mat_solve_block_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b, slong block_size, flint_rand_t state)
{
    int good = 0, ret;
    mp_ptr x1;
    nmod_sparse_vec_t z;
    nmod_sparse_mat_t Mb;
    if (M->r != M->c) return 0; /* TODO */
    if (_nmod_vec_is_zero(b, M->c))
    {
        _nmod_vec_zero(x, M->c);
        return 1;
    }

    /* TODO: Precondition M */
    x1 = _nmod_vec_init(M->c + 1);
    nmod_sparse_vec_init(z);
    nmod_sparse_mat_init(Mb, M->r, M->c, M->mod);
    nmod_sparse_mat_set(Mb, M);
    nmod_sparse_mat_append_col(Mb, b);
    nmod_sparse_mat_append_row(Mb, z);

    ret = nmod_sparse_mat_nullvector_block_wiedemann(x1, Mb, block_size, state);
    if (ret && x1[M->c] != UWORD(0)) 
    {
        _nmod_vec_scalar_mul_nmod(x, x1, M->c, nmod_neg(nmod_inv(x1[M->c], M->mod), M->mod), M->mod);
        good = 1;
    }
    nmod_sparse_vec_clear(z);
    nmod_sparse_mat_clear(Mb);
    _nmod_vec_clear(x1);
    return good;
}

int nmod_sparse_mat_nullvector_block_wiedemann(mp_ptr x, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state) 
{
    int ret = 0;
    slong l, ns, k;
    slong *d;
    mp_ptr b;
    nmod_mat_struct Y[3], *S;
    if (M->r != M->c) return 0; /* TODO */

    ns = 2*M->r/block_size + 3; /* Maybe 5? */
    S = flint_malloc(ns*sizeof(*S));
    d = flint_calloc(2*block_size, sizeof(*d));
    b = _nmod_vec_init(M->r);
    for (k = 0; k < ns; ++k) nmod_mat_init(&S[k], block_size, block_size, M->mod.n);
    for (l = 0; l < 3; ++l) nmod_mat_init(&Y[l], M->c, block_size, M->mod.n);
    do nmod_mat_randfull(&Y[0], state);
    while (nmod_mat_is_zero(&Y[0]));

    nmod_sparse_mat_mul_mat(&Y[1], M, &Y[0]);
    make_block_sequences(S, ns, M, &Y[1]);
    find_block_min_poly(S, d, ns, M->r);

    for (l = 0; l < block_size; ++l)
    {
        nmod_mat_set(&Y[1], &Y[0]);
        make_block_sum(x, S, d, M, Y + 1, l);
        nmod_sparse_mat_mul_vec(b, M, x);
        if (!_nmod_vec_is_zero(x, M->c) && _nmod_vec_is_zero(b, M->r)) {ret = 1; break;};
    }    
    return ret;
}


