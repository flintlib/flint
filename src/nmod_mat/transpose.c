/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "machine_vectors.h"
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

static void
_transpose_blocked(nn_ptr B, slong Bstride, nn_srcptr A, slong Astride, slong m, slong n)
{
    slong i, j;

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

    if (m <= 8)     /* few rows of A: each line of A serves 8 rows of B */
    {
        for (i = 0; i < n; i++)
            for (j = 0; j < m; j++)
                B[i * Bstride + j] = A[j * Astride + i];
        return;
    }

    if (n <= 8)     /* few columns of A: each line of B serves 8 rows of A */
    {
        for (j = 0; j < m; j++)
            for (i = 0; i < n; i++)
                B[i * Bstride + j] = A[j * Astride + i];
        return;
    }

    if (m <= TRANSPOSE_BLOCK && n <= TRANSPOSE_BLOCK)
    {
        _transpose_basecase(B, Bstride, A, Astride, m, n);
        return;
    }

    for (j = 0; j < m; j += TRANSPOSE_BLOCK)
        for (i = 0; i < n; i += TRANSPOSE_BLOCK)
            _transpose_basecase(B + i * Bstride + j, Bstride, A + j * Astride + i, Astride,
                                FLINT_MIN(TRANSPOSE_BLOCK, m - j),
                                FLINT_MIN(TRANSPOSE_BLOCK, n - i));
}

#if defined(VEC4N_TRANSPOSE)

/*
    With vector instructions (AVX2, NEON), 4 x 4 blocks are transposed in
    registers by VEC4N_TRANSPOSE. Out of place, they are grouped in 8 x 8
    blocks, so that each row of A and of B that a block touches is 64
    contiguous bytes, and the blocks are swept along strips of 8 rows of B,
    restricted to 64 rows of A at a time: when the rows are not 64-byte
    aligned, the lines of A that a strip uses only in part are then still in
    L1 for the next strip. With fewer than 8 rows or columns, the 4 x 4
    blocks are used alone. In place, pairs of 4 x 4 blocks are swapped. What
    remains (fewer than 8, or 4, rows or columns) is done entrywise.
*/
#define TRANSPOSE_SIMD_ROWS 64

/*
    Write rows of B by 32-byte vector stores slows things down on Intel
    (measured on Ice Lake and Sapphire Rapids) up to 2x,
    -> two 16-byte stores avoid this
    But on AMD 32-byte stores are similar or faster (mesured on zen 4)...
*/
#if defined(__znver3__) || defined(__znver4__) || defined(__znver5__)
# define _transpose_store vec4n_store_unaligned
#else
# define _transpose_store vec4n_store_unaligned_halves
#endif

/* B = transpose of A, for 4 x 4 blocks; A and B may be equal */
FLINT_FORCE_INLINE void
_transpose_4x4(nn_ptr B, slong Bstride, nn_srcptr A, slong Astride)
{
    vec4n a0 = vec4n_load_unaligned(A);
    vec4n a1 = vec4n_load_unaligned(A + Astride);
    vec4n a2 = vec4n_load_unaligned(A + 2 * Astride);
    vec4n a3 = vec4n_load_unaligned(A + 3 * Astride);

    VEC4N_TRANSPOSE(a0, a1, a2, a3, a0, a1, a2, a3);

    _transpose_store(B, a0);
    _transpose_store(B + Bstride, a1);
    _transpose_store(B + 2 * Bstride, a2);
    _transpose_store(B + 3 * Bstride, a3);
}

/* X, Y = transpose of Y, transpose of X, for 4 x 4 blocks of stride stride */
FLINT_FORCE_INLINE void
_transpose_swap_4x4(nn_ptr X, nn_ptr Y, slong stride)
{
    vec4n x0 = vec4n_load_unaligned(X);
    vec4n x1 = vec4n_load_unaligned(X + stride);
    vec4n x2 = vec4n_load_unaligned(X + 2 * stride);
    vec4n x3 = vec4n_load_unaligned(X + 3 * stride);
    vec4n y0 = vec4n_load_unaligned(Y);
    vec4n y1 = vec4n_load_unaligned(Y + stride);
    vec4n y2 = vec4n_load_unaligned(Y + 2 * stride);
    vec4n y3 = vec4n_load_unaligned(Y + 3 * stride);

    VEC4N_TRANSPOSE(x0, x1, x2, x3, x0, x1, x2, x3);
    VEC4N_TRANSPOSE(y0, y1, y2, y3, y0, y1, y2, y3);

    _transpose_store(Y, x0);
    _transpose_store(Y + stride, x1);
    _transpose_store(Y + 2 * stride, x2);
    _transpose_store(Y + 3 * stride, x3);
    _transpose_store(X, y0);
    _transpose_store(X + stride, y1);
    _transpose_store(X + 2 * stride, y2);
    _transpose_store(X + 3 * stride, y3);
}

#endif

void
_nmod_mat_transpose(nn_ptr B, slong Bstride, nn_srcptr A, slong Astride, slong m, slong n)
{
#if defined(VEC4N_TRANSPOSE)
    slong i, j, j0, j1, k, m8 = m - m % 8, n8 = n - n % 8;

    if (m < 8 || n < 8)     /* 4 x 4 blocks */
    {
        slong m4 = m - m % 4, n4 = n - n % 4;

        for (i = 0; i < n4; i += 4)
            for (j = 0; j < m4; j += 4)
                _transpose_4x4(B + i * Bstride + j, Bstride, A + j * Astride + i, Astride);

        if (n > n4)
            _transpose_blocked(B + n4 * Bstride, Bstride, A + n4, Astride, m, n - n4);
        if (m > m4)
            _transpose_blocked(B + m4, Bstride, A + m4 * Astride, Astride, m - m4, n4);
        return;
    }

    for (j0 = 0; j0 < m8; j0 = j1)
    {
        j1 = FLINT_MIN(j0 + TRANSPOSE_SIMD_ROWS, m8);

        for (i = 0; i < n8; i += 8)
        {
            for (j = j0; j < j1; j += 8)
            {
                _transpose_4x4(B + i * Bstride + j, Bstride, A + j * Astride + i, Astride);
                _transpose_4x4(B + i * Bstride + j + 4, Bstride, A + (j + 4) * Astride + i, Astride);
                _transpose_4x4(B + (i + 4) * Bstride + j, Bstride, A + j * Astride + i + 4, Astride);
                _transpose_4x4(B + (i + 4) * Bstride + j + 4, Bstride, A + (j + 4) * Astride + i + 4, Astride);
            }

            if (j1 == m8)   /* the last m % 8 rows of A, on these columns */
                for (j = m8; j < m; j++)
                    for (k = i; k < i + 8; k++)
                        B[k * Bstride + j] = A[j * Astride + k];
        }

        /* the last n % 8 columns of A, on these rows */
        if (n > n8)
            _transpose_blocked(B + n8 * Bstride + j0, Bstride, A + j0 * Astride + n8, Astride, j1 - j0, n - n8);
    }

    if (m > m8 && n > n8)
        _transpose_blocked(B + n8 * Bstride + m8, Bstride, A + m8 * Astride + n8, Astride, m - m8, n - n8);
#else
    _transpose_blocked(B, Bstride, A, Astride, m, n);
#endif
}

/* A (n x n) = transpose of A */
static void
_nmod_mat_transpose_inplace(nn_ptr A, slong stride, slong n)
{
    slong i0, j0, i, j;

#if defined(VEC4N_TRANSPOSE)
    /* below 8 x 8, the loop below is as fast or faster (2x at 4 x 4 on M4) */
    if (n >= 8)
    {
        slong k, n4 = n - n % 4;

        for (i = 0; i < n4; i += 4)
        {
            _transpose_4x4(A + i * stride + i, stride, A + i * stride + i, stride);

            for (j = i + 4; j < n4; j += 4)
                _transpose_swap_4x4(A + i * stride + j, A + j * stride + i, stride);

            for (k = i; k < i + 4; k++)     /* the last n % 4 columns */
                for (j = n4; j < n; j++)
                    FLINT_SWAP(ulong, A[k * stride + j], A[j * stride + k]);
        }

        for (k = n4; k < n; k++)            /* the bottom-right corner */
            for (j = k + 1; j < n; j++)
                FLINT_SWAP(ulong, A[k * stride + j], A[j * stride + k]);

        return;
    }
#endif

    /* swap the tiles (i0, j0) and (j0, i0) */
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
