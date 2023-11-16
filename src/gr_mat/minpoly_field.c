/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"
#include "gr_mat.h"

int gr_mat_reduce_row(slong * column, gr_mat_t A, slong * P, slong * L, slong m, gr_ctx_t ctx);

int
gr_mat_minpoly_field(gr_poly_t p, const gr_mat_t X, gr_ctx_t ctx)
{
    slong n = X->r, i, j, c, c1, c2, r1, r2;
    slong  * P1, * P2, * L1, * L2;
    gr_mat_t A, B, v;
    int first_poly = 1, indep = 1;
    gr_poly_t b, g, r;
    gr_ptr t, h;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    TMP_INIT;

    if (X->r != X->c)
        return GR_DOMAIN;

    if (n == 0)
        return gr_poly_one(p, ctx);

    if (n == 1)
    {
        gr_poly_fit_length(p, 2, ctx);
        status |= gr_neg(GR_ENTRY(p->coeffs, 0, sz), GR_MAT_ENTRY(X, 0, 0, sz), ctx);
        status |= gr_one(GR_ENTRY(p->coeffs, 1, sz), ctx);
        _gr_poly_set_length(p, 2, ctx);
        return status;
    }

    TMP_START;

    GR_TMP_INIT2(t, h, ctx);

    gr_init(h, ctx);
    gr_poly_init(b, ctx);
    gr_poly_init(g, ctx);
    gr_poly_init(r, ctx);
    gr_mat_init(A, n + 1, 2*n + 1, ctx);
    gr_mat_init(B, n, n, ctx);
    gr_mat_init(v, n, 1, ctx);

    status |= gr_poly_one(p, ctx);

    L1 = (slong *) TMP_ALLOC((n + 1)*sizeof(slong));
    L2 = (slong *) TMP_ALLOC(n*sizeof(slong));
    P1 = (slong *) TMP_ALLOC((2*n + 1)*sizeof(slong));
    P2 = (slong *) TMP_ALLOC(n*sizeof(slong));

    for (i = 1; i <= n + 1; i++)
        L1[i - 1] = n + i;

    for (i = 1; i <= n; i++)
        L2[i - 1] = n;

    for (i = 1; i < n; i++)
        P2[i] = -1;

    P2[0] = 0;

    r2 = c2 = 0;
    first_poly = 1;

    while (r2 < n)
    {
        for (i = 0; i < 2*n + 1; i++)
            P1[i] = -1;

        for (i = 0; i < n; i++)
        {
            status |= gr_zero(GR_MAT_ENTRY(v, i, 0, sz), ctx);
            status |= gr_zero(GR_MAT_ENTRY(B, r2, i, sz), ctx);
            status |= gr_zero(GR_MAT_ENTRY(A, 0, i, sz), ctx);
        }

        P1[c2] = 0;
        P2[c2] = r2;

        status |= gr_one(GR_MAT_ENTRY(v, c2, 0, sz), ctx);
        status |= gr_one(GR_MAT_ENTRY(B, r2, c2, sz), ctx);
        status |= gr_one(GR_MAT_ENTRY(A, 0, c2, sz), ctx);
        status |= gr_one(GR_MAT_ENTRY(A, 0, n, sz), ctx);

        indep = 1;

        r1 = 0;
        c1 = -1;

        while (c1 < n && r1 < n)
        {
            r1++;
            r2 = indep ? r2 + 1 : r2;

            status |= gr_mat_mul(v, X, v, ctx);

            for (i = 0; i < n; i++)
                status |= gr_set(GR_MAT_ENTRY(A, r1, i, sz), GR_MAT_ENTRY(v, i, 0, sz), ctx);

            for (i = n; i < n + r1; i++)
                status |= gr_zero(GR_MAT_ENTRY(A, r1, i, sz), ctx);

            status |= gr_one(GR_MAT_ENTRY(A, r1, n + r1, sz), ctx);
            if (status != GR_SUCCESS)
                goto cleanup;

            status |= gr_mat_reduce_row(&c1, A, P1, L1, r1, ctx);
            if (status != GR_SUCCESS)
                goto cleanup;

            if (indep && r2 < n && !first_poly)
            {
                for (i = 0; i < n; i++)
                    status |= gr_set(GR_MAT_ENTRY(B, r2, i, sz), GR_MAT_ENTRY(v, i, 0, sz), ctx);

                status |= gr_mat_reduce_row(&c, B, P2, L2, r2, ctx);
                if (status != GR_SUCCESS)
                    goto cleanup;

                indep = (c != -1);
            }
        }

        if (first_poly)
        {
            for (i = 0; i < n; i++)
                P2[i] = P1[i];

            r2 = r1;
        }

        c = -1;

        for (i = c2 + 1; i < n; i++)
        {
            if (P2[i] == -1)
            {
                c = i;
                break;
            }
        }

        c2 = c;

        gr_poly_fit_length(b, r1 + 1, ctx);
        status |= gr_inv(h, GR_MAT_ENTRY(A, r1, n + r1, sz), ctx);
        for (i = 0; i < r1 + 1; i++)
        {
            status |= gr_mul(t, GR_MAT_ENTRY(A, r1, n + i, sz), h, ctx);
            status |= gr_poly_set_coeff_scalar(b, i, t, ctx);
        }
        _gr_poly_set_length(b, r1 + 1, ctx);

        status |= gr_poly_gcd(g, p, b, ctx);
        status |= gr_poly_mul(p, p, b, ctx);
        status |= gr_poly_divrem(p, r, p, g, ctx);

        if (first_poly && r2 < n)
        {
            for (i = 0; i < r1; i++)
            {
                for (j = 0; j < n; j++)
                    status |= gr_set(GR_MAT_ENTRY(B, i, j, sz),  GR_MAT_ENTRY(A, i, j, sz), ctx);
            }
        }

        first_poly = 0;
    }

cleanup:
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(v, ctx);

    gr_poly_clear(b, ctx);
    gr_poly_clear(g, ctx);
    gr_poly_clear(r, ctx);

    GR_TMP_CLEAR2(t, h, ctx);

    TMP_END;

    return status;
}
