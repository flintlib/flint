/*
    Copyright (C) 2010,2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

/* gcc11 on x86_64 generates significantly worse code with -O3 */
#pragma GCC optimize("-O2")

/*
with op = 0, computes D = A*B
with op = 1, computes D = C + A*B
with op = -1, computes D = C - A*B
*/

static inline void
_nmod_mat_addmul_basic_op(nn_ptr * D, nn_ptr * const C, nn_ptr * const A,
    nn_ptr * const B, slong m, slong k, slong n, int op, nmod_t mod, dot_params_t params)
{
    slong i, j;
    ulong c;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            c = _nmod_vec_dot_ptr(A[i], B, j, k, mod, params);

            if (op == 1)
                c = nmod_add(C[i][j], c, mod);
            else if (op == -1)
                c = nmod_sub(C[i][j], c, mod);

            D[i][j] = c;
        }
    }
}

static inline void
_nmod_mat_addmul_transpose_op(nn_ptr * D, const nn_ptr * C, const nn_ptr * A,
    const nn_ptr * B, slong m, slong k, slong n, int op, nmod_t mod, dot_params_t params)
{
    nn_ptr tmp;
    ulong c;
    slong i, j;

    tmp = flint_malloc(sizeof(ulong) * k * n);

    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
            tmp[j*k + i] = B[i][j];

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            c = _nmod_vec_dot(A[i], tmp + j*k, k, mod, params);

            if (op == 1)
                c = nmod_add(C[i][j], c, mod);
            else if (op == -1)
                c = nmod_sub(C[i][j], c, mod);

            D[i][j] = c;
        }
    }

    flint_free(tmp);
}

/* Assumes nlimbs = 1 */
static void
_nmod_mat_addmul_packed_op(nn_ptr * D, const nn_ptr * C, const nn_ptr * A,
    const nn_ptr * B, slong M, slong N, slong K, int op, nmod_t mod)
{
    slong i, j, k;
    slong Kpack;
    int pack, pack_bits;
    ulong c, d, mask;
    nn_ptr tmp;
    nn_ptr Aptr, Tptr;

    /* bound unreduced entry */
    c = N * (mod.n-1) * (mod.n-1);
    pack_bits = FLINT_BIT_COUNT(c);
    pack = FLINT_BITS / pack_bits;
    Kpack = (K + pack - 1) / pack;

    if (pack_bits == FLINT_BITS)
        mask = UWORD(-1);
    else
        mask = (UWORD(1) << pack_bits) - 1;

    tmp = _nmod_vec_init(Kpack * N);

    /* pack and transpose B */
    for (i = 0; i < Kpack; i++)
    {
        for (k = 0; k < N; k++)
        {
            c = B[k][i * pack];

            for (j = 1; j < pack && i * pack + j < K; j++)
                c |= B[k][i * pack + j] << (pack_bits * j);

            tmp[i * N + k] = c;
        }
    }

    /* multiply */
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < Kpack; j++)
        {
            Aptr = A[i];
            Tptr = tmp + j * N;

            c = 0;

            /* unroll by 4 */
            for (k = 0; k + 4 <= N; k += 4)
            {
                c += Aptr[k + 0] * Tptr[k + 0];
                c += Aptr[k + 1] * Tptr[k + 1];
                c += Aptr[k + 2] * Tptr[k + 2];
                c += Aptr[k + 3] * Tptr[k + 3];
            }

            for ( ; k < N; k++)
                c += Aptr[k] * Tptr[k];

            /* unpack and reduce */
            for (k = 0; k < pack && j * pack + k < K; k++)
            {
                d = (c >> (k * pack_bits)) & mask;
                NMOD_RED(d, d, mod);

                if (op == 1)
                    d = nmod_add(C[i][j * pack + k], d, mod);
                else if (op == -1)
                    d = nmod_sub(C[i][j * pack + k], d, mod);

                D[i][j * pack + k] = d;
            }
        }
    }

    _nmod_vec_clear(tmp);
}



void
_nmod_mat_mul_classical_op(nmod_mat_t D, const nmod_mat_t C,
                                const nmod_mat_t A, const nmod_mat_t B, int op)
{
    slong m, k, n;
    nmod_t mod;

    mod = A->mod;
    m = A->r;
    k = A->c;
    n = B->c;

    if (k == 0 || mod.n == 1)  // covers params.method == _DOT0
    {
        if (op == 0)
            nmod_mat_zero(D);
        else
            nmod_mat_set(D, C);
        return;
    }

    const dot_params_t params = _nmod_vec_dot_params(k, mod);

    // TODO vec_dot changes --> thresholds to re-examine
    if (params.method == _DOT1 && m > 10 && k > 10 && n > 10)
    {
        _nmod_mat_addmul_packed_op(D->rows, (op == 0) ? NULL : C->rows,
            A->rows, B->rows, m, k, n, op, D->mod);
    }
    else if (m < NMOD_MAT_MUL_TRANSPOSE_CUTOFF
        || n < NMOD_MAT_MUL_TRANSPOSE_CUTOFF
        || k < NMOD_MAT_MUL_TRANSPOSE_CUTOFF)
    {
        _nmod_mat_addmul_basic_op(D->rows, (op == 0) ? NULL : C->rows,
            A->rows, B->rows, m, k, n, op, D->mod, params);
    }
    else
    {
        _nmod_mat_addmul_transpose_op(D->rows, (op == 0) ? NULL : C->rows,
            A->rows, B->rows, m, k, n, op, D->mod, params);
    }
}


void
nmod_mat_mul_classical(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    _nmod_mat_mul_classical_op(C, NULL, A, B, 0);
}
