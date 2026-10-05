/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"

#define TRANSPOSE_BLOCK 32

/* entrywise: B (n x m) = transpose of A (m x n) */
static void
_transpose_basecase(nn_ptr B, slong Bstride, nn_srcptr A, slong Astride, slong m, slong n)
{
    slong i, j;

    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
            B[i * Bstride + j] = A[j * Astride + i];
}

void
_nmod_mat_transpose(nn_ptr B, slong Bstride, nn_srcptr A, slong Astride, slong m, slong n)
{
    slong i, j;

    if (m <= TRANSPOSE_BLOCK && n <= TRANSPOSE_BLOCK)
    {
        _transpose_basecase(B, Bstride, A, Astride, m, n);
        return;
    }

    if (m == 1)     /* row to column */
    {
        for (i = 0; i < n; i++)
            B[i * Bstride] = A[i];
        return;
    }

    if (n == 1)     /* column to row */
    {
        for (j = 0; j < m; j++)
            B[j] = A[j * Astride];
        return;
    }

    for (j = 0; j < m; j += TRANSPOSE_BLOCK)
        for (i = 0; i < n; i += TRANSPOSE_BLOCK)
            _transpose_basecase(B + i * Bstride + j, Bstride, A + j * Astride + i, Astride,
                                FLINT_MIN(TRANSPOSE_BLOCK, m - j),
                                FLINT_MIN(TRANSPOSE_BLOCK, n - i));
}

/* A (n x n) = transpose of A, swapping the tiles (i0, j0) and (j0, i0) */
static void
_nmod_mat_transpose_inplace(nn_ptr A, slong stride, slong n)
{
    slong i0, j0, i, j;

    for (i0 = 0; i0 < n; i0 += TRANSPOSE_BLOCK)
        for (j0 = i0; j0 < n; j0 += TRANSPOSE_BLOCK)
            for (i = i0; i < FLINT_MIN(i0 + TRANSPOSE_BLOCK, n); i++)
                for (j = FLINT_MAX(j0, i + 1); j < FLINT_MIN(j0 + TRANSPOSE_BLOCK, n); j++)
                    FLINT_SWAP(ulong, A[i * stride + j], A[j * stride + i]);
}

void
nmod_mat_transpose(nmod_mat_t B, const nmod_mat_t A)
{
    if (B->r != A->c || B->c != A->r)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_mat_transpose). Incompatible dimensions.\n");
    }

    if (A == B) /* must be square */
        _nmod_mat_transpose_inplace(A->entries, A->stride, A->r);
    else
        _nmod_mat_transpose(B->entries, B->stride, A->entries, A->stride, A->r, A->c);
}
