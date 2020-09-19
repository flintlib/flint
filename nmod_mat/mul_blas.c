/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 William Hart
    This file is part of FLINT.
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"
#include "thread_support.h"
#include "profiler.h"

#if HAVE_BLAS

#include "cblas.h"

/*
    b (and hence the signed remainder) definitely fits into an int, which
    helps the conversion to a double
*/
static void _lift_vec(double * a, ulong * b, slong len, ulong n)
{
    slong i;
    for (i = 0; i < len; i++)
        a[i] = (int)(b[i] - (n & FLINT_SIGN_EXT(n/2 - b[i])));
}

/*
    the shift for negative a must satisfy
        ctx.n divides shift, and
        2^53 <= shift <= 2^64
*/
static mp_limb_t _reduce(slong a, nmod_t ctx, ulong shift)
{
    mp_limb_t b, c = (a < 0) ? a + shift : a;
    NMOD_RED(b, c, ctx);
    return b;
}

int nmod_mat_mul_blas(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, j;
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    double * dC, * dA, * dB;
    ulong hi, lo, shift;
    nmod_t ctx = C->mod;
    timeit_t timer;

    if (m < 1 || k < 1 || n < 1 || m > INT_MAX || k > INT_MAX || n > INT_MAX)
        return 0;

    /*
        Each nmod has |_lift()| <= floor(mod.n/2) in the signed representation.
        Want k*floor(mod.n/2)^2 < 2^53 **assuming** no BLAS rounding errors.
    */

    if (FLINT_BITS != 64)
        return 0;

    umul_ppmm(hi, lo, ctx.n/2, ctx.n/2);
    if (hi != 0)
        return 0;

    umul_ppmm(hi, lo, lo, k);
    if (hi != 0 || lo >= (UWORD(1) << 53))
        return 0;

    dA = flint_malloc(m*k*sizeof(double));
    dB = flint_malloc(k*n*sizeof(double));
    dC = flint_calloc(m*n, sizeof(double));

    for (i = 0; i < m; i++)
        _lift_vec(dA + i*k, A->rows[i], k, ctx.n);

    for (i = 0; i < k; i++)
        _lift_vec(dB + i*n, B->rows[i], n, ctx.n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                                                1.0, dA, k, dB, n, 0.0, dC, n);

    shift = ((UWORD(1)<<54)/ctx.n)*ctx.n;
    for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
        C->rows[i][j] = _reduce(dC[i*n + j], ctx, shift);

    flint_free(dA);
    flint_free(dB);
    flint_free(dC);

    return 1;
}

#else

int nmod_mat_mul_blas(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    return 0;
}

#endif
