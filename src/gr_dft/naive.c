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
#include "gr_dft.h"

/* res must not alias vec. Scrambled order for the naive algorithm is
   bit-reversed, matching GR_DFT_ALG_CT. */
int
_gr_dft_naive(gr_ptr res, gr_srcptr vec, int inverse, int scrambled,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong n = P->n, j, k, e;
    gr_ptr t, rtmp = NULL;

    GR_TMP_INIT(t, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    for (k = 0; k < n; k++)
    {
        gr_ptr acc = GR_ENTRY(res,
                (!inverse && scrambled) ? n_revbin(k, P->depth) : k, sz);

        e = 0;
        status |= gr_set(acc,
                GR_ENTRY(vec, (inverse && scrambled) ? n_revbin(0, P->depth) : 0, sz), ctx);

        for (j = 1; j < n; j++)
        {
            gr_srcptr src = GR_ENTRY(vec,
                    (inverse && scrambled) ? n_revbin(j, P->depth) : j, sz);

            e += k;
            if (e >= n)
                e -= n;
            status |= _gr_dft_mul_root(t, src, e, inverse, rtmp, P);
            status |= gr_add(acc, acc, t, ctx);
        }
    }

    GR_TMP_CLEAR(t, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    return status;
}
