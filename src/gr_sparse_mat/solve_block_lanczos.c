/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 Kartik Venkatram

    Algorithm taken from P. Montgomery, "A Block Lanczos Algorithm for 
    Finding Dependencies over GF(2)", Advances in Cryptology - EUROCRYPT '95

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

/* Run row Gaussian elimination on first b block of [T I], save that, */
/* if no pivot is found for a given column c, we Gaussian eliminate   */
/* the column c + b and zero out the row c. In addition, we reorder   */
/* columns so that ones corresponding to zero entries in S go first.  */
/* See Figure 1 in the above reference for details. */
static int compute_nWi_S(slong *rk, gr_mat_t nWi, int *S, const gr_mat_t Torig, gr_ctx_t ctx) 
{
    const slong b = Torig->r;
    slong pc, i, j;
    gr_mat_t T;
    gr_mat_struct *X;
    slong *P;
    gr_ptr cc;
    int status = GR_SUCCESS;

    GR_TMP_INIT(cc, ctx);
    P = flint_malloc(b * sizeof(*P));
    gr_mat_init(T, b, b, ctx);
    status |= gr_mat_set(T, Torig, ctx);
    status |= gr_mat_one(nWi, ctx);

    /* Set permutation to have previously dependent vectors at front */
    P = flint_malloc(b*sizeof(*P));
    j = 0;
    for (i = 0; i < b; ++i) if (!S[i]) P[j++] = i;
    for (i = 0; i < b; ++i) if (S[i]) P[j++] = i;
    
    *rk = 0;
    for (j = 0; j < b; ++j) 
    {
        pc = P[j]; /* Pivot col */

        /* Find viable pivot row (from T if possible, then from W) */
        for (X = T, i = j; i < b; ++i)
            if (gr_is_zero(gr_mat_entry_srcptr(X, P[i], pc, ctx), ctx) == T_FALSE)
                break;
        if (i == b)
            for (X = nWi, i = j; i < b; ++i)
                if (gr_is_zero(gr_mat_entry_srcptr(X, P[i], pc, ctx), ctx) == T_FALSE)
                    break;
        S[pc] = X == T; /* Viable column in V */
        status |= gr_mat_swap_rows(T, NULL, pc, P[i], ctx);
        status |= gr_mat_swap_rows(nWi, NULL, pc, P[i], ctx); /* Now pivot row = pivot col */

        /* Make pivot one */
        status |= gr_inv(cc, gr_mat_entry_ptr(X, pc, pc, ctx), ctx);
        status |= _gr_vec_mul_scalar(T->rows[pc], T->rows[pc], b, cc, ctx);
        status |= _gr_vec_mul_scalar(nWi->rows[pc], nWi->rows[pc], b, cc, ctx);


        /* Kill all other entries in pivot column */
        for (i = 0; i < b; ++i)
        {
            status |= gr_neg(cc, gr_mat_entry_ptr(X, P[i], pc, ctx), ctx);
            if (i == j || cc == 0) continue;
            status |= _gr_vec_addmul_scalar(T->rows[P[i]], T->rows[pc], T->c, cc, ctx);
            status |= _gr_vec_addmul_scalar(nWi->rows[P[i]], nWi->rows[pc], nWi->c, cc, ctx);
        }
        if (S[pc]) (*rk)++; /* Count viable columns */
        else 
        {
            /* Kill row of both matrices */
            status |= _gr_vec_zero(T->rows[pc], b, ctx);
            status |= _gr_vec_zero(nWi->rows[pc], b, ctx);
        }
    }

    status |= gr_mat_neg(nWi, nWi, ctx);
    gr_mat_clear(T, ctx);
    GR_TMP_CLEAR(cc, ctx);

    return status;
}

static int kill_columns(gr_mat_t M, int *good, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j;
    for (j = 0; j < M->c; ++j)
        if (good[j] == 0)
            for (i = 0; i < M->r; ++i)
                status |= gr_zero(gr_mat_entry_ptr(M, i, j, ctx), ctx);
    return status;
}

int gr_lil_mat_solve_block_lanczos(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, slong block_size, flint_rand_t state, gr_ctx_t ctx) 
{
    slong r, c, i, j, prev_i, next_i, iter, cur_dim, total_dim = 0;
    gr_lil_mat_t Mt; /* Transpose of M, we work with A = MtM */
    gr_mat_struct V[3]; /* Keep track of current vector and two previous ones */
    gr_mat_t MV; /* Application of M to V */
    gr_mat_t AV; /* Application of Mt to MV */
    int *SSt; /* S is the maximal projection s.t. (VS)^tAVS is invertible, so SSt kills the dropped columns */
    gr_mat_struct nWi[3]; /* -S((VS)^tAVS)^-1S^t */
    gr_mat_t VSSt; /* V with invalid vectors zeroed out */
    gr_mat_t T; /* Used to store transposes for inner products */
    gr_mat_t VtAV; /* Inner product <V, V>_A */
    gr_mat_t AVtAVSSt_VtAV; /* Sum <AV, V>_A SS^t + <V, V>_A, shared by two updates */
    gr_mat_t DEF; /* Used to store coefficient matrices D, E, and F */
    gr_mat_t I, tmp; /* I_{b x b}, tmp used as scratch */
    gr_ptr Mtb, SStVtMtb, WiSStVtMtb, VSStWiSStVtMtb; /* Intermediate elements in x update */
    int status = GR_SUCCESS;

    // TODO: handle this
    if (x == b)
        return GR_DOMAIN;

    r = gr_sparse_mat_nrows(M, ctx);
    c = gr_sparse_mat_ncols(M, ctx);

    if (_gr_vec_is_zero(b, r, ctx))
    {
        status |= _gr_vec_zero(x, c, ctx);
        return GR_SUCCESS;
    }

    gr_lil_mat_init(Mt, c, r, ctx);
    for (i = 0; i < 3; ++i)
        gr_mat_init(&V[i], c, block_size, ctx);
    gr_mat_init(MV, r, block_size, ctx); /* Intermediate product */
    gr_mat_init(AV, c, block_size, ctx); /* Symmetric product */
    SSt = flint_malloc(block_size*sizeof(*SSt));
    for (i = 0; i < 3; ++i)
        gr_mat_init(&nWi[i], block_size, block_size, ctx);
    gr_mat_init(VSSt, c, block_size, ctx);
    gr_mat_init(T, block_size, c, ctx); /* Transpose for computing matrix dot products */
    gr_mat_init(VtAV, block_size, block_size, ctx);
    gr_mat_init(AVtAVSSt_VtAV, block_size, block_size, ctx); /* (AV)^T(AV) + VtAV */
    gr_mat_init(DEF, block_size, block_size, ctx); /* Shared by D, E, and F */
    gr_mat_init(I, block_size, block_size, ctx);
    gr_mat_init(tmp, block_size, block_size, ctx);
    GR_TMP_INIT_VEC(Mtb, c, ctx);
    GR_TMP_INIT_VEC(SStVtMtb, block_size, ctx);
    GR_TMP_INIT_VEC(WiSStVtMtb, block_size, ctx);
    GR_TMP_INIT_VEC(VSStWiSStVtMtb, c, ctx);

    status |= _gr_vec_zero(x, c, ctx);
    status |= gr_lil_mat_transpose(Mt, M, ctx);
    for (i = 0; i < block_size; ++i) SSt[i] = 1;
    status |= gr_mat_one(I, ctx);
    status |= gr_lil_mat_mul_vec(Mtb, Mt, b, ctx);

    /* Initialize V[0] randomly */
    for (i = 0; i < c; ++i)
        for (j = 0; j < block_size; ++j)
            status |= gr_randtest(gr_mat_entry_ptr(&V[0], i, j, ctx), state, ctx);

    for (iter = 0; ; ++iter) 
    {
        i = iter % 3;
        next_i = (iter + 1) % 3;
        prev_i = (iter + 2) % 3;
        if (iter >= 2)
        {
            /* Compute the F value for this round (minus the final term) */
            status |= gr_mat_addmul(DEF, I, VtAV, &nWi[prev_i], ctx);
            status |= gr_mat_mul(tmp, &nWi[next_i], DEF, ctx);
            status |= gr_mat_mul(DEF, tmp, AVtAVSSt_VtAV, ctx);
        }

        /* Compute AV and V'AV */
        status |= gr_lil_mat_mul_mat(MV, M, &V[i], ctx);
        status |= gr_lil_mat_mul_mat(AV, Mt, MV, ctx);
        status |= gr_mat_transpose(T, &V[i], ctx);
        status |= gr_mat_mul(VtAV, T, AV, ctx);
        if (gr_mat_is_zero(VtAV, ctx) == T_TRUE) {status = GR_SUCCESS; break;}

        /* Compute W^{-1} and indices of bad vectors */
        status |= compute_nWi_S(&cur_dim, &nWi[i], SSt, VtAV, ctx);
        total_dim += cur_dim;
        if (cur_dim == 0 || total_dim > c) break; /* Ran out of vectors */

        /* Update x_i = x_{i-1} - (VSS^t) W^{-1} (VSS^t)^tb */
        status |= gr_mat_set(VSSt, &V[i], ctx);
        status |= kill_columns(VSSt, SSt, ctx);
        status |= gr_mat_transpose(T, VSSt, ctx);
        status |= gr_mat_mul_vec(SStVtMtb, T, Mtb, ctx);
        status |= gr_mat_mul_vec(WiSStVtMtb, &nWi[i], SStVtMtb, ctx);
        status |= gr_mat_mul_vec(VSStWiSStVtMtb, VSSt, WiSStVtMtb, ctx);
        status |= _gr_vec_add(x, x, VSStWiSStVtMtb, c, ctx);

        /**
         * Per Equation (19), we compute the next vector 
         *   V_{i+1} = AV_iS_iS_i^t + V_i D + V_{i-1} E + V_{i-2} F 
         * where
         *   D = I - W_i^-1((AV_i)^tAV_iS_iS_i^t + V_i^tAV_i)
         *   E = -W_{i-1}^-1V_i^tAV_iS_iS_i^t
         *   F = -W_{i-2}^-1(I - V_{i-1}^tAV_{i-1}W_{i-1}^-1)
         *                  ((AV_{i-1})^tAV_{i-1}S_{i-1}S_{i-1}^t + V_{i-1}^tAV_{i-1})S_iS_i^t
         **/
        if (iter >= 2)
        {
            /* V_{i+1} = V_{i-2} F */
            status |= kill_columns(DEF, SSt, ctx);
            status |= gr_mat_mul(VSSt, &V[next_i], DEF, ctx);
            status |= gr_mat_set(&V[next_i], VSSt, ctx);
        }
        if (iter >= 1)
        {
            /* V_{i+1} += V_{i-1} E */
            status |= gr_mat_mul(DEF, &nWi[prev_i], VtAV, ctx);
            status |= kill_columns(DEF, SSt, ctx);
            status |= gr_mat_addmul(&V[next_i], &V[next_i], &V[prev_i], DEF, ctx);
        }
        /* V_{i+1} += V_i D */
        status |= gr_mat_transpose(T, AV, ctx);
        status |= gr_mat_mul(tmp, T, AV, ctx);
        status |= kill_columns(tmp, SSt, ctx);
        status |= gr_mat_add(AVtAVSSt_VtAV, tmp, VtAV, ctx);
        status |= gr_mat_addmul(DEF, I, &nWi[i], AVtAVSSt_VtAV, ctx);
        status |= gr_mat_addmul(&V[next_i], &V[next_i], &V[i], DEF, ctx);

        /* V_{i+1} += AVSS^t */
        status |= kill_columns(AV, SSt, ctx);
        status |= gr_mat_add(&V[next_i], &V[next_i], AV, ctx);

        if (gr_mat_is_zero(&V[next_i], ctx) == T_TRUE) {status = GR_SUCCESS; break;}
    }
    status |= _gr_vec_neg(x, x, c, ctx);
    gr_lil_mat_clear(Mt, ctx);
    for (i = 0; i < 3; ++i)
        gr_mat_clear(&V[i], ctx);
    gr_mat_clear(MV, ctx);
    gr_mat_clear(AV, ctx);
    flint_free(SSt);
    for (i = 0; i < 3; ++i)
        gr_mat_clear(&nWi[i], ctx);
    gr_mat_clear(T, ctx);
    gr_mat_clear(VtAV, ctx);
    gr_mat_clear(VSSt, ctx);
    gr_mat_clear(AVtAVSSt_VtAV, ctx);
    gr_mat_clear(DEF, ctx);
    gr_mat_clear(I, ctx);
    gr_mat_clear(tmp, ctx);
    GR_TMP_CLEAR_VEC(Mtb, c, ctx);
    GR_TMP_CLEAR_VEC(SStVtMtb, block_size, ctx);
    GR_TMP_CLEAR_VEC(WiSStVtMtb, block_size, ctx);
    GR_TMP_CLEAR_VEC(VSStWiSStVtMtb, c, ctx);
    return status;
}

int gr_lil_mat_nullvector_block_lanczos(gr_ptr x, const gr_lil_mat_t M, slong block_size, flint_rand_t state, gr_ctx_t ctx) 
{
    int status = GR_SUCCESS;
    gr_ptr x2, b;
    GR_TMP_INIT_VEC(x2, M->c, ctx);
    GR_TMP_INIT_VEC(b, M->r, ctx);

    status |= _gr_vec_randtest(x, state, M->c, ctx);
    status |= gr_lil_mat_mul_vec(b, M, x, ctx);
    status |= gr_lil_mat_solve_block_lanczos(x2, M, b, block_size, state, ctx);

    if (status == GR_SUCCESS) 
    {
        status |= _gr_vec_sub(x, x, x2, M->c, ctx);
        if (_gr_vec_is_zero(x, M->c, ctx) == T_TRUE)
            status = GR_TEST_FAIL;
        else
        {
            status |= gr_lil_mat_mul_vec(b, M, x, ctx);
            if (_gr_vec_is_zero(b, M->r, ctx) == T_FALSE)
                status = GR_DOMAIN;
        }
    }
    GR_TMP_CLEAR_VEC(x2, M->c, ctx);
    GR_TMP_CLEAR_VEC(b, M->r, ctx);
    return status;
}
