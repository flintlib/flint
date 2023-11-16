/*
    Copyright (C) 2013 Timo Hartmann
    Copyright (C) 2018, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

/* todo: optimize; use dot/norm functions */
/* todo: numerical improvements */
int
gr_mat_hessenberg_householder(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong n;
    gr_ptr T;
    gr_ptr H, G, F, f, ff;
    slong i, j, k;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    n = mat->r;
    if (n != mat->c)
        return GR_DOMAIN;

    status |= gr_mat_set(res, mat, ctx);

#define MAT_ENTRY(i, j) GR_MAT_ENTRY(res, i, j, sz)

    GR_TMP_INIT_VEC(T, n, ctx);
    GR_TMP_INIT5(H, G, F, f, ff, ctx);

    for (i = n - 1; i >= 2; i--)
    {
        status |= gr_zero(H, ctx);

        for (k = 0; k < i; k++)
        {
            status |= gr_conj(ff, MAT_ENTRY(i, k), ctx);
            status |= gr_addmul(H, MAT_ENTRY(i, k), ff, ctx);
        }

        status |= gr_set(F, MAT_ENTRY(i, i - 1), ctx);
        status |= gr_abs(f, F, ctx);
        status |= gr_sqrt(G, H, ctx);
        status |= gr_neg(MAT_ENTRY(i, i - 1), G, ctx);
        status |= gr_div(ff, F, f, ctx);
        status |= gr_mul(GR_ENTRY(T, i, sz), G, ff, ctx);
        status |= gr_add(GR_ENTRY(T, i, sz), GR_ENTRY(T, i, sz), F, ctx);
        status |= gr_mul(MAT_ENTRY(i, i - 1), MAT_ENTRY(i, i - 1), ff, ctx);
        status |= gr_addmul(H, G, f, ctx);
        status |= gr_rsqrt(H, H, ctx);
        status |= gr_mul(GR_ENTRY(T, i, sz), GR_ENTRY(T, i, sz), H, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        for (k = 0; k < i - 1; k++)
            status |= gr_mul(MAT_ENTRY(i, k), MAT_ENTRY(i, k), H, ctx);

        for (j = 0; j < i; j++)
        {
            status |= gr_conj(ff, GR_ENTRY(T, i, sz), ctx);
            status |= gr_mul(G, ff, MAT_ENTRY(j, i - 1), ctx);

            for (k = 0; k < i - 1; k++)
            {
                status |= gr_conj(ff, MAT_ENTRY(i, k), ctx);
                status |= gr_addmul(G, ff, MAT_ENTRY(j, k), ctx);
            }

            status |= gr_submul(MAT_ENTRY(j, i - 1), G, GR_ENTRY(T, i, sz), ctx);
            for (k = 0; k < i - 1; k++)
                status |= gr_submul(MAT_ENTRY(j, k), G, MAT_ENTRY(i, k), ctx);
        }

        for (j = 0; j < n; j++)
        {
            status |= gr_mul(G, GR_ENTRY(T, i, sz), MAT_ENTRY(i - 1, j), ctx);
            for (k = 0; k < i - 1; k++)
                status |= gr_addmul(G, MAT_ENTRY(i, k), MAT_ENTRY(k, j), ctx);

            status |= gr_conj(ff, GR_ENTRY(T, i, sz), ctx);
            status |= gr_submul(MAT_ENTRY(i - 1, j), G, ff, ctx);
            for (k = 0; k < i - 1; k++)
            {
                status |= gr_conj(ff, MAT_ENTRY(i, k), ctx);
                status |= gr_submul(MAT_ENTRY(k, j), G, ff, ctx);
            }
        }
    }

    for (i = 0; i < n; i++)
        for (j = i + 2; j < n; j++)
            status |= gr_zero(MAT_ENTRY(j, i), ctx);

cleanup:
    GR_TMP_CLEAR_VEC(T, n, ctx);
    GR_TMP_CLEAR5(H, G, F, f, ff, ctx);

    return status;
}
