/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"

/* Bluestein's algorithm for odd n. With t = (n+1)/2 = 1/2 mod n we have
   j k = t (j^2 + k^2 - (j-k)^2) mod n, so

       X_k = w^(t k^2) sum_j (w^(t j^2) x_j) w^(-t (k-j)^2),

   a cyclic convolution of length n of the chirped input with the even,
   n-periodic kernel h_j = w^(-t j^2), embedded in a cyclic convolution
   of power-of-two length conv_len >= 2n - 1 which is evaluated with
   the sub-plan P->P1. The chirp factors are entries of the root table,
   so in complex Karatsuba mode they enjoy the reduced-multiplication classes.
   The transformed kernel P->bl_kern is precomputed with the 1/conv_len
   scaling of the inverse sub-transform folded in, and is stored in the
   (scrambled) output order of the sub-plan, so no reordering passes
   are needed.

   The inverse transform is the forward transform applied to the
   cyclically reversed input, x'_j = x_((n - j) mod n).

   res must not alias vec. */
int
_gr_dft_bluestein(gr_ptr res, gr_srcptr vec, int inverse,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong rsz = (P->real_ctx != NULL) ? P->real_ctx->sizeof_elem : 0;
    ulong n = P->n, N = P->conv_len, j, k, e, t2, ninv;
    gr_ptr bufA, bufB, rtmp = NULL;

    bufA = flint_malloc(2 * N * sz);
    _gr_vec_init(bufA, 2 * N, ctx);
    bufB = GR_ENTRY(bufA, N, sz);

    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    t2 = (n + 1) / 2;
    ninv = n_preinvert_limb(n);

    /* chirp the (possibly reversed) input */
    for (j = 0; j < n; j++)
    {
        gr_srcptr src = GR_ENTRY(vec,
                (inverse && j > 0) ? n - j : j, sz);

        e = n_mulmod2_preinv(j, j, n, ninv);
        e = n_mulmod2_preinv(e, t2, n, ninv);
        status |= _gr_dft_mul_root(GR_ENTRY(bufA, j, sz), src, e, 0, rtmp, P);
    }

    status |= _gr_vec_zero(GR_ENTRY(bufA, n, sz), N - n, ctx);

    /* convolve with the kernel */
    status |= _gr_dft_precomp_raw(bufB, bufA, 0, P->P1, ctx);

    for (k = 0; k < N; k++)
        status |= _gr_dft_mul_const(GR_ENTRY(bufA, k, sz),
                GR_ENTRY(bufB, k, sz),
                GR_ENTRY(P->bl_kern, k, sz),
                (P->bl_wtab != NULL) ? GR_ENTRY(P->bl_wtab, 3 * k, rsz) : NULL,
                rtmp, P);

    status |= _gr_dft_precomp_raw(bufB, bufA, 1, P->P1, ctx);

    /* chirp the output */
    for (k = 0; k < n; k++)
    {
        e = n_mulmod2_preinv(k, k, n, ninv);
        e = n_mulmod2_preinv(e, t2, n, ninv);
        status |= _gr_dft_mul_root(GR_ENTRY(res, k, sz),
                GR_ENTRY(bufB, k, sz), e, 0, rtmp, P);
    }

    /* over the fixed-point contexts the kernel carries an extra
       factor 1/2 (see the kernel construction); undo it with exact
       doublings */
    if (P->bl_shifted)
        for (k = 0; k < n; k++)
            status |= gr_add(GR_ENTRY(res, k, sz), GR_ENTRY(res, k, sz),
                    GR_ENTRY(res, k, sz), ctx);

    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    _gr_vec_clear(bufA, 2 * N, ctx);
    flint_free(bufA);

    return status;
}
