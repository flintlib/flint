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

#ifdef T

#include <string.h>
#include "templates.h"

/* Run row Gaussian elimination on first b block of [T I], save that, */
/* if no pivot is found for a given column c, we Gaussian eliminate   */
/* the column c + b and zero out the row c. In addition, we reorder   */
/* columns so that ones corresponding to zero entries in S go first.  */
/* See Figure 1 in the above reference for details. */
static int compute_nWi_S(TEMPLATE(T, mat_t) nWi, int *S, const TEMPLATE(T, mat_t) Torig, const TEMPLATE(T, ctx_t) ctx) 
{
    
    const slong b = Torig->r;
    slong pc, i, j, rk = 0;
    slong *P;
    TEMPLATE(T, t) cc;
    TEMPLATE(T, mat_t) T;
    TEMPLATE(T, mat_struct) *X;

    P = flint_malloc(b * sizeof(*P));
    TEMPLATE(T, init) (cc, ctx);
    TEMPLATE(T, mat_init) (T, b, b, ctx);
    TEMPLATE(T, mat_set) (T, Torig, ctx);
    TEMPLATE(T, mat_one) (nWi, ctx);

    /* Set permutation to have previously dependent vectors at front */
    P = flint_malloc(b*sizeof(*P));
    j = 0;
    for (i = 0; i < b; ++i) if (!S[i]) P[j++] = i;
    for (i = 0; i < b; ++i) if (S[i]) P[j++] = i;
    
    for (j = 0; j < b; ++j) 
    {
        pc = P[j]; /* Pivot col */

        /* Find viable pivot row (from T if possible, then from W) */
        for (X = T, i = j; i < b && TEMPLATE(T, is_zero)(&X->rows[P[i]][pc], ctx); ++i); 
        if (i == b)
            for (X = nWi, i = j; i < b && TEMPLATE(T, is_zero)(&X->rows[P[i]][pc], ctx); ++i); 
        S[pc] = X == T; /* Viable column in V */
        TEMPLATE(T, mat_swap_rows) (T, NULL, pc, P[i], ctx);
        TEMPLATE(T, mat_swap_rows) (nWi, NULL, pc, P[i], ctx); /* Now pivot row = pivot col */

        /* Make pivot one */
        TEMPLATE(T, inv) (cc, &X->rows[pc][pc], ctx);
        _TEMPLATE(T, TEMPLATE(vec_scalar_mul, T)) (T->rows[pc], T->rows[pc], b, cc, ctx);
        _TEMPLATE(T, TEMPLATE(vec_scalar_mul, T)) (nWi->rows[pc], nWi->rows[pc], b, cc, ctx);

        /* Kill all other entries in pivot column */
        for (i = 0; i < b; ++i)
        {
            TEMPLATE(T, neg) (cc, &X->rows[P[i]][pc], ctx);
            if (i == j || TEMPLATE(T, is_zero) (cc, ctx)) continue;
            _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T)) (T->rows[P[i]], T->rows[pc], T->c, cc, ctx);
            _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T)) (nWi->rows[P[i]], nWi->rows[pc], nWi->c, cc, ctx);
        }
        if (S[pc]) rk++; /* Count viable columns */
        else 
        {
            /* Kill row of both matrices */
            _TEMPLATE(T, vec_zero) (T->rows[pc], b, ctx);
            _TEMPLATE(T, vec_zero) (nWi->rows[pc], b, ctx);
        }
    }
    TEMPLATE(T, mat_neg) (nWi, nWi, ctx);
    TEMPLATE(T, mat_clear) (T, ctx);

    return rk;
}

static void kill_columns(TEMPLATE(T, mat_t) M, int *good, const TEMPLATE(T, ctx_t) ctx)
{
    slong r, c;
    for (c = 0; c < M->c; ++c)
        if (good[c] == 0)
            for (r = 0; r < M->r; ++r)
                TEMPLATE(T, zero) (&M->rows[r][c], ctx);
}

int TEMPLATE(T, sparse_mat_solve_block_lanczos) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *b, slong block_size, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx) 
{
    int ret = 0;
    slong i, prev_i, next_i, iter, cur_dim, total_dim = 0;
    TEMPLATE(T, sparse_mat_t) Mt; /* Transpose of M, we work with A = MtM */
    TEMPLATE(T, mat_struct) V[3]; /* Keep track of current vector and two previous ones */
    TEMPLATE(T, mat_t) MV; /* Application of M to V */
    TEMPLATE(T, mat_t) AV; /* Application of Mt to MV */
    int *SSt; /* S is the maximal projection s.t. (VS)^tAVS is invertible, so SSt kills the dropped columns */
    TEMPLATE(T, mat_struct) nWi[3]; /* -S((VS)^tAVS)^-1S^t */
    TEMPLATE(T, mat_t) VSSt; /* V with invalid vectors zeroed out */
    TEMPLATE(T, mat_t) T; /* Used to store transposes for inner products */
    TEMPLATE(T, mat_t) VtAV; /* Inner product <V, V>_A */
    TEMPLATE(T, mat_t) AVtAVSSt_VtAV; /* Sum <AV, V>_A SS^t + <V, V>_A, shared by two updates */
    TEMPLATE(T, mat_t) DEF; /* Used to store coefficient matrices D, E, and F */
    TEMPLATE(T, mat_t) I, tmp; /* I_{b x b} */
    TEMPLATE(T, struct) *Mtb, *SStVtMtb, *WiSStVtMtb, *VSStWiSStVtMtb; /* Intermediate elements in x update */
    if (_TEMPLATE(T, vec_is_zero) (b, M->r, ctx))
    {
        _TEMPLATE(T, vec_zero) (x, M->c, ctx);
        return 1;
    }
    TEMPLATE(T, sparse_mat_init) (Mt, M->c, M->r, ctx);
    for (i = 0; i < 3; ++i) TEMPLATE(T, mat_init) (&V[i], M->c, block_size, ctx);
    TEMPLATE(T, mat_init) (MV, M->r, block_size, ctx); /* Intermediate product */
    TEMPLATE(T, mat_init) (AV, M->c, block_size, ctx); /* Symmetric product */
    SSt = flint_malloc(block_size*sizeof(*SSt));
    for (i = 0; i < 3; ++i) TEMPLATE(T, mat_init) (&nWi[i], block_size, block_size, ctx);
    TEMPLATE(T, mat_init) (VSSt, M->c, block_size, ctx); /* Transpose for computing matrix dot products */
    TEMPLATE(T, mat_init) (T, block_size, M->c, ctx); /* Transpose for computing matrix dot products */
    TEMPLATE(T, mat_init) (VtAV, block_size, block_size, ctx);
    TEMPLATE(T, mat_init) (AVtAVSSt_VtAV, block_size, block_size, ctx); /* (AV)^T(AV) + VtAV */
    TEMPLATE(T, mat_init) (DEF, block_size, block_size, ctx); /* Shared by D, E, and F */
    TEMPLATE(T, mat_init) (I, block_size, block_size, ctx);
    TEMPLATE(T, mat_init) (tmp, block_size, block_size, ctx);
    Mtb = _TEMPLATE(T, vec_init) (M->c, ctx);
    SStVtMtb = _TEMPLATE(T, vec_init) (block_size, ctx);
    WiSStVtMtb = _TEMPLATE(T, vec_init) (block_size, ctx);
    VSStWiSStVtMtb = _TEMPLATE(T, vec_init) (M->c, ctx);

    _TEMPLATE(T, vec_zero) (x, M->c, ctx);
    TEMPLATE(T, sparse_mat_transpose) (Mt, M, ctx);
    for (i = 0; i < block_size; ++i) SSt[i] = 1;
    TEMPLATE(T, mat_one) (I, ctx);
    TEMPLATE(T, sparse_mat_mul_vec) (Mtb, Mt, b, ctx);

    /* Initialize V[0] randomly */
    for (i = 0; i < V[0].r*V[0].c; ++i)
        TEMPLATE(T, randtest) (&V[0].entries[i], state, ctx);

    for (iter = 0; ; ++iter) 
    {
        i = iter % 3;
        next_i = (iter + 1) % 3;
        prev_i = (iter + 2) % 3;
        if (iter >= 2)
        {
            /* Compute the F value for this round (minus the final term) */
            TEMPLATE(T, mat_addmul) (DEF, I, VtAV, &nWi[prev_i], ctx);
            TEMPLATE(T, mat_mul) (tmp, &nWi[next_i], DEF, ctx);
            TEMPLATE(T, mat_mul) (DEF, tmp, AVtAVSSt_VtAV, ctx);
        }

        /* Compute AV and V'AV */
        TEMPLATE(T, sparse_mat_mul_mat) (MV, M, &V[i], ctx);
        TEMPLATE(T, sparse_mat_mul_mat) (AV, Mt, MV, ctx);
        TEMPLATE(T, mat_transpose) (T, &V[i], ctx);
        TEMPLATE(T, mat_mul) (VtAV, T, AV, ctx);
        if (TEMPLATE(T, mat_is_zero) (VtAV, ctx)) {ret = 1; break;}
        
        /* Compute W^{-1} and indices of bad vectors */
        cur_dim = compute_nWi_S(&nWi[i], SSt, VtAV, ctx);
        total_dim += cur_dim;
        if (cur_dim == 0 || total_dim > M->c) break; /* Ran out of vectors */

        /* Update x_i = x_{i-1} - (VSS^t) W^{-1} (VSS^t)^tb */
        TEMPLATE(T, mat_set) (VSSt, &V[i], ctx);
        kill_columns(VSSt, SSt, ctx);
        TEMPLATE(T, mat_transpose) (T, VSSt, ctx);
        TEMPLATE(T, mat_mul_vec) (SStVtMtb, T, Mtb, ctx);
        TEMPLATE(T, mat_mul_vec) (WiSStVtMtb, &nWi[i], SStVtMtb, ctx);
        TEMPLATE(T, mat_mul_vec) (VSStWiSStVtMtb, VSSt, WiSStVtMtb, ctx);
        _TEMPLATE(T, vec_add) (x, x, VSStWiSStVtMtb, M->c, ctx);

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
            kill_columns(DEF, SSt, ctx);
            TEMPLATE(T, mat_mul) (VSSt, &V[next_i], DEF, ctx);
            TEMPLATE(T, mat_set) (&V[next_i], VSSt, ctx);
        }
        if (iter >= 1)
        {
            /* V_{i+1} += V_{i-1} E */
            TEMPLATE(T, mat_mul) (DEF, &nWi[prev_i], VtAV, ctx);
            kill_columns(DEF, SSt, ctx);
            TEMPLATE(T, mat_addmul) (&V[next_i], &V[next_i], &V[prev_i], DEF, ctx);
        }
        /* V_{i+1} += V_i D */
        TEMPLATE(T, mat_transpose) (T, AV, ctx);
        TEMPLATE(T, mat_mul) (tmp, T, AV, ctx);
        kill_columns(tmp, SSt, ctx);
        TEMPLATE(T, mat_add) (AVtAVSSt_VtAV, tmp, VtAV, ctx);
        TEMPLATE(T, mat_addmul) (DEF, I, &nWi[i], AVtAVSSt_VtAV, ctx);
        TEMPLATE(T, mat_addmul) (&V[next_i], &V[next_i], &V[i], DEF, ctx);

        /* V_{i+1} += AVSS^t */
        kill_columns(AV, SSt, ctx);
        TEMPLATE(T, mat_add) (&V[next_i], &V[next_i], AV, ctx);

        if (TEMPLATE(T, mat_is_zero) (&V[i], ctx)) {ret = 1; break;}
    }
    _TEMPLATE(T, vec_neg) (x, x, M->c, ctx);
    TEMPLATE(T, sparse_mat_clear) (Mt, ctx);
    for (i = 0; i < 3; ++i) TEMPLATE(T, mat_clear) (&V[i], ctx);
    TEMPLATE(T, mat_clear) (MV, ctx);
    TEMPLATE(T, mat_clear) (AV, ctx);
    flint_free(SSt);
    for (i = 0; i < 3; ++i) TEMPLATE(T, mat_clear) (&nWi[i], ctx);
    TEMPLATE(T, mat_clear) (T, ctx);
    TEMPLATE(T, mat_clear) (VtAV, ctx);
    TEMPLATE(T, mat_clear) (VSSt, ctx);
    TEMPLATE(T, mat_clear) (AVtAVSSt_VtAV, ctx);
    TEMPLATE(T, mat_clear) (DEF, ctx);
    TEMPLATE(T, mat_clear) (I, ctx);
    TEMPLATE(T, mat_clear) (tmp, ctx);
    _TEMPLATE(T, vec_clear) (Mtb, M->c, ctx);
    _TEMPLATE(T, vec_clear) (SStVtMtb, block_size, ctx);
    _TEMPLATE(T, vec_clear) (WiSStVtMtb, block_size, ctx);
    _TEMPLATE(T, vec_clear) (VSStWiSStVtMtb, M->c, ctx);
    return ret;
}

int TEMPLATE(T, sparse_mat_nullvector_block_lanczos) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, slong block_size, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx) 
{
    int ret = 1;
    TEMPLATE(T, struct) *x2, *b;
    x2 = _TEMPLATE(T, vec_init) (M->c, ctx);
    b = _TEMPLATE(T, vec_init) (M->r, ctx);

    _TEMPLATE(T, vec_randtest) (x, state, M->c, ctx);
    TEMPLATE(T, sparse_mat_mul_vec) (b, M, x, ctx);
    if (TEMPLATE(T, sparse_mat_solve_block_lanczos) (x2, M, b, block_size, state, ctx) == 0) ret = 0; /* Lanczos failed */
    if (ret)
    {
        _TEMPLATE(T, vec_sub) (x, x, x2, M->c, ctx);
        TEMPLATE(T, sparse_mat_mul_vec) (b, M, x, ctx);
        ret = !_TEMPLATE(T, vec_is_zero) (x, M->c, ctx) && _TEMPLATE(T, vec_is_zero) (b, M->r, ctx);
    }
    _TEMPLATE(T, vec_clear) (x2, M->c, ctx);
    _TEMPLATE(T, vec_clear) (b, M->r, ctx);
    return ret;
}


#endif
