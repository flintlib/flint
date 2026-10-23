/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"

/* Good-Thomas prime factor algorithm for n = n1 n2 with gcd(n1, n2) = 1.
   With the input index map j = (n2 j1 + n1 j2) mod n,

       w^(j k) = (w^n2)^(j1 k1) (w^n1)^(j2 k2),

   where k1 = k mod n1 and k2 = k mod n2, so the transform is a pure
   two-dimensional DFT (no twiddle factors) with roots w^n2 (rows of
   length n1, sub-plan P1) and w^n1 (columns of length n2, sub-plan P2),
   followed by the CRT output map k = (pfa_a k1 + pfa_b k2) mod n. All
   index maps are evaluated incrementally.

   res must not alias vec. */
int
_gr_dft_pfa(gr_ptr res, gr_srcptr vec, int inverse,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong n = P->n, n1 = P->n1, n2 = P->n2;
    ulong j1, j2, k1, k2, jb, j, kb, k;
    ulong L = FLINT_MAX(n1, n2);
    gr_ptr mat, buf, buf2;

    mat = flint_malloc((n + 2 * L) * sz);
    _gr_vec_init(mat, n + 2 * L, ctx);
    buf = GR_ENTRY(mat, n, sz);
    buf2 = GR_ENTRY(buf, L, sz);

    /* gather: mat[j1 n2 + j2] = vec[(n2 j1 + n1 j2) mod n] */
    for (j1 = 0, jb = 0; j1 < n1; j1++)
    {
        for (j2 = 0, j = jb; j2 < n2; j2++)
        {
            status |= gr_set(GR_ENTRY(mat, j1 * n2 + j2, sz),
                    GR_ENTRY(vec, j, sz), ctx);
            j += n1;
            if (j >= n)
                j -= n;
        }
        jb += n2;
        if (jb >= n)
            jb -= n;
    }

    /* rows: n1 transforms of length n2 with sub-plan P2 */
    for (j1 = 0; j1 < n1; j1++)
    {
        gr_ptr row = GR_ENTRY(mat, j1 * n2, sz);
        status |= _gr_dft_precomp_raw(buf, row, inverse, P->P2, ctx);
        status |= _gr_vec_set(row, buf, n2, ctx);
    }

    /* columns: n2 transforms of length n1 with sub-plan P1 */
    for (j2 = 0; j2 < n2; j2++)
    {
        for (j1 = 0; j1 < n1; j1++)
            status |= gr_set(GR_ENTRY(buf, j1, sz),
                    GR_ENTRY(mat, j1 * n2 + j2, sz), ctx);

        status |= _gr_dft_precomp_raw(buf2, buf, inverse, P->P1, ctx);

        for (j1 = 0; j1 < n1; j1++)
            status |= gr_set(GR_ENTRY(mat, j1 * n2 + j2, sz),
                    GR_ENTRY(buf2, j1, sz), ctx);
    }

    /* scatter: res[(pfa_a k1 + pfa_b k2) mod n] = mat[k1 n2 + k2] */
    for (k1 = 0, kb = 0; k1 < n1; k1++)
    {
        for (k2 = 0, k = kb; k2 < n2; k2++)
        {
            status |= gr_set(GR_ENTRY(res, k, sz),
                    GR_ENTRY(mat, k1 * n2 + k2, sz), ctx);
            k += P->pfa_b;
            if (k >= n)
                k -= n;
        }
        kb += P->pfa_a;
        if (kb >= n)
            kb -= n;
    }

    _gr_vec_clear(mat, n + 2 * L, ctx);
    flint_free(mat);

    return status;
}
