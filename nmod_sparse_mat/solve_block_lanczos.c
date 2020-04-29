/*
    Copyright (C) 2020 Kartik Venkatram

    Algorithm taken from P. Montgomery, "A Block Lanczos Algorithm for 
    Finding Dependencies over GF(2)", Advances in Cryptology - EUROCRYPT '95

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

/* Run row Gaussian elimination on first b block of [T I], save that, */
/* if no pivot is found for a given column c, we Gaussian eliminate   */
/* the column c + b and zero out the row c. In addition, we reorder   */
/* columns so that ones corresponding to zero entries in S go first.  */
/* See Figure 1 in the above reference for details. */
static int compute_nWi_S(nmod_mat_t nWi, int *S, const nmod_mat_t Torig) 
{
    
    const slong b = Torig->r;
    slong pc, i, j, rk = 0;
    nmod_mat_t T;
    nmod_mat_struct *X;
    slong *P;
    mp_limb_t cc;

    P = flint_malloc(b * sizeof(*P));
    nmod_mat_init(T, b, b, Torig->mod.n);
    nmod_mat_set(T, Torig);
    nmod_mat_one(nWi);

    /* Set permutation to have previously dependent vectors at front */
    P = flint_malloc(b*sizeof(*P));
    j = 0;
    for (i = 0; i < b; ++i) if (!S[i]) P[j++] = i;
    for (i = 0; i < b; ++i) if (S[i]) P[j++] = i;
    
    for (j = 0; j < b; ++j) 
    {
        pc = P[j]; /* Pivot col */

        /* Find viable pivot row (from T if possible, then from W) */
        for (X = T, i = j; i < b && X->rows[P[i]][pc] == 0; ++i);
        if (i == b)
            for (X = nWi, i = j; i < b && X->rows[P[i]][pc] == 0; ++i);
        S[pc] = X == T; /* Viable column in V */
        nmod_mat_swap_rows(T, NULL, pc, P[i]);
        nmod_mat_swap_rows(nWi, NULL, pc, P[i]); /* Now pivot row = pivot col */

        /* Make pivot one */
        cc = nmod_inv(X->rows[pc][pc], T->mod);
        _nmod_vec_scalar_mul_nmod(T->rows[pc], T->rows[pc], b, cc, T->mod);
        _nmod_vec_scalar_mul_nmod(nWi->rows[pc], nWi->rows[pc], b, cc, T->mod);


        /* Kill all other entries in pivot column */
        for (i = 0; i < b; ++i)
        {
            cc = nmod_neg(X->rows[P[i]][pc], T->mod);
            if (i == j || cc == 0) continue;
            _nmod_vec_scalar_addmul_nmod(T->rows[P[i]], T->rows[pc], T->c, cc, T->mod);
            _nmod_vec_scalar_addmul_nmod(nWi->rows[P[i]], nWi->rows[pc], nWi->c, cc, T->mod);
        }
        if (S[pc]) rk++; /* Count viable columns */
        else 
        {
            /* Kill row of both matrices */
            _nmod_vec_zero(T->rows[pc], b);
            _nmod_vec_zero(nWi->rows[pc], b);
        }
    }

    nmod_mat_neg(nWi, nWi);
    nmod_mat_clear(T);

    return rk;
}

static void kill_columns(nmod_mat_t M, int *good)
{
    slong r, c;
    for (c = 0; c < M->c; ++c)
        if (good[c] == 0)
            for (r = 0; r < M->r; ++r)
                M->rows[r][c] = UWORD(0);
}

int nmod_sparse_mat_solve_block_lanczos(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b, slong block_size, flint_rand_t state) 
{
    int ret = 0;
    slong i, prev_i, next_i, iter, cur_dim, total_dim = 0;
    nmod_sparse_mat_t Mt; /* Transpose of M, we work with A = MtM */
    nmod_mat_struct V[3]; /* Keep track of current vector and two previous ones */
    nmod_mat_t MV; /* Application of M to V */
    nmod_mat_t AV; /* Application of Mt to MV */
    int *SSt; /* S is the maximal projection s.t. (VS)^tAVS is invertible, so SSt kills the dropped columns */
    nmod_mat_struct nWi[3]; /* -S((VS)^tAVS)^-1S^t */
    nmod_mat_t VSSt; /* V with invalid vectors zeroed out */
    nmod_mat_t T; /* Used to store transposes for inner products */
    nmod_mat_t VtAV; /* Inner product <V, V>_A */
    nmod_mat_t AVtAVSSt_VtAV; /* Sum <AV, V>_A SS^t + <V, V>_A, shared by two updates */
    nmod_mat_t DEF; /* Used to store coefficient matrices D, E, and F */
    nmod_mat_t I, tmp; /* I_{b x b}, tmp used as scratch */
    mp_ptr Mtb, SStVtMtb, WiSStVtMtb, VSStWiSStVtMtb; /* Intermediate elements in x update */

    if (_nmod_vec_is_zero(b, M->r))
    {
        _nmod_vec_zero(x, M->c);
        return 1;
    }

    nmod_sparse_mat_init(Mt, M->c, M->r, M->mod);
    for (i = 0; i < 3; ++i) nmod_mat_init(&V[i], M->c, block_size, M->mod.n);
    nmod_mat_init(MV, M->r, block_size, M->mod.n); /* Intermediate product */
    nmod_mat_init(AV, M->c, block_size, M->mod.n); /* Symmetric product */
    SSt = flint_malloc(block_size*sizeof(*SSt));
    for (i = 0; i < 3; ++i) nmod_mat_init(&nWi[i], block_size, block_size, M->mod.n);
    nmod_mat_init(VSSt, M->c, block_size, M->mod.n);
    nmod_mat_init(T, block_size, M->c, M->mod.n); /* Transpose for computing matrix dot products */
    nmod_mat_init(VtAV, block_size, block_size, M->mod.n);
    nmod_mat_init(AVtAVSSt_VtAV, block_size, block_size, M->mod.n); /* (AV)^T(AV) + VtAV */
    nmod_mat_init(DEF, block_size, block_size, M->mod.n); /* Shared by D, E, and F */
    nmod_mat_init(I, block_size, block_size, M->mod.n);
    nmod_mat_init(tmp, block_size, block_size, M->mod.n);
    Mtb = _nmod_vec_init(M->c);
    SStVtMtb = _nmod_vec_init(block_size);
    WiSStVtMtb = _nmod_vec_init(block_size);
    VSStWiSStVtMtb = _nmod_vec_init(M->c);

    _nmod_vec_zero(x, M->c);
    nmod_sparse_mat_transpose(Mt, M);
    for (i = 0; i < block_size; ++i) SSt[i] = 1;
    nmod_mat_one(I);
    nmod_sparse_mat_mul_vec(Mtb, Mt, b);

    /* Initialize V[0] randomly */
    for (i = 0; i < V[0].r*V[0].c; ++i)
        V[0].entries[i] = n_randint(state, V[0].mod.n);

    for (iter = 0; ; ++iter) 
    {
        i = iter % 3;
        next_i = (iter + 1) % 3;
        prev_i = (iter + 2) % 3;
        if (iter >= 2)
        {
            /* Compute the F value for this round (minus the final term) */
            nmod_mat_addmul(DEF, I, VtAV, &nWi[prev_i]);
            nmod_mat_mul(tmp, &nWi[next_i], DEF);
            nmod_mat_mul(DEF, tmp, AVtAVSSt_VtAV);
        }

        /* Compute AV and V'AV */
        nmod_sparse_mat_mul_mat(MV, M, &V[i]);
        nmod_sparse_mat_mul_mat(AV, Mt, MV);
        nmod_mat_transpose(T, &V[i]);
        nmod_mat_mul(VtAV, T, AV);
        if (nmod_mat_is_zero(VtAV)) {ret = 1; break;}

        /* Compute W^{-1} and indices of bad vectors */
        cur_dim = compute_nWi_S(&nWi[i], SSt, VtAV);
        total_dim += cur_dim;
        if (cur_dim == 0 || total_dim > M->c) break; /* Ran out of vectors */

        /* Update x_i = x_{i-1} - (VSS^t) W^{-1} (VSS^t)^tb */
        nmod_mat_set(VSSt, &V[i]);
        kill_columns(VSSt, SSt);
        nmod_mat_transpose(T, VSSt);
        nmod_mat_mul_vec(SStVtMtb, T, Mtb);
        nmod_mat_mul_vec(WiSStVtMtb, &nWi[i], SStVtMtb);
        nmod_mat_mul_vec(VSStWiSStVtMtb, VSSt, WiSStVtMtb);
        _nmod_vec_add(x, x, VSStWiSStVtMtb, M->c, M->mod);

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
            kill_columns(DEF, SSt);
            nmod_mat_mul(VSSt, &V[next_i], DEF);
            nmod_mat_set(&V[next_i], VSSt);
        }
        if (iter >= 1)
        {
            /* V_{i+1} += V_{i-1} E */
            nmod_mat_mul(DEF, &nWi[prev_i], VtAV);
            kill_columns(DEF, SSt);
            nmod_mat_addmul(&V[next_i], &V[next_i], &V[prev_i], DEF);
        }
        /* V_{i+1} += V_i D */
        nmod_mat_transpose(T, AV);
        nmod_mat_mul(tmp, T, AV);
        kill_columns(tmp, SSt);
        nmod_mat_add(AVtAVSSt_VtAV, tmp, VtAV);
        nmod_mat_addmul(DEF, I, &nWi[i], AVtAVSSt_VtAV);
        nmod_mat_addmul(&V[next_i], &V[next_i], &V[i], DEF);

        /* V_{i+1} += AVSS^t */
        kill_columns(AV, SSt);
        nmod_mat_add(&V[next_i], &V[next_i], AV);

        if (nmod_mat_is_zero(&V[next_i])) {ret = 1; break;}
    }
    _nmod_vec_neg(x, x, M->c, M->mod);
    nmod_sparse_mat_clear(Mt);
    for (i = 0; i < 3; ++i) nmod_mat_clear(&V[i]);
    nmod_mat_clear(MV);
    nmod_mat_clear(AV);
    flint_free(SSt);
    for (i = 0; i < 3; ++i) nmod_mat_clear(&nWi[i]);
    nmod_mat_clear(T);
    nmod_mat_clear(VtAV);
    nmod_mat_clear(VSSt);
    nmod_mat_clear(AVtAVSSt_VtAV);
    nmod_mat_clear(DEF);
    nmod_mat_clear(I);
    nmod_mat_clear(tmp);
    _nmod_vec_clear(SStVtMtb);
    _nmod_vec_clear(WiSStVtMtb);
    _nmod_vec_clear(VSStWiSStVtMtb);
    _nmod_vec_clear(Mtb);
    return ret;
}

int nmod_sparse_mat_nullvector_block_lanczos(mp_ptr x, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state) 
{
    int ret = 1;
    mp_ptr x2, b;
    x2 = _nmod_vec_init(M->c);
    b = _nmod_vec_init(M->r);

    _nmod_vec_randtest(x, state, M->c, M->mod);
    nmod_sparse_mat_mul_vec(b, M, x);
    if (nmod_sparse_mat_solve_block_lanczos(x2, M, b, block_size, state) == 0) ret = 0; /* Lanczos failed */
    if (ret)
    {
        _nmod_vec_sub(x, x, x2, M->c, M->mod);
        nmod_sparse_mat_mul_vec(b, M, x);
        ret = !_nmod_vec_is_zero(x, M->c) && _nmod_vec_is_zero(b, M->r);
    }
    _nmod_vec_clear(x2);
    _nmod_vec_clear(b);
    return ret;
}