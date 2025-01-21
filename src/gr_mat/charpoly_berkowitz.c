/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"

int
_gr_mat_charpoly_berkowitz(gr_ptr cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    slong n = mat->r;
    int status = GR_SUCCESS;

#define MAT(ii, jj) GR_MAT_ENTRY(mat, ii, jj, sz)
#define POL(ii) GR_ENTRY(cp, ii, sz)

    if (n == 0)
    {
        status |= gr_one(POL(0), ctx);
    }
    else if (n == 1)
    {
        status |= gr_neg(POL(0), MAT(0, 0), ctx);
        status |= gr_one(POL(1), ctx);
    }
    else if (n == 2)
    {
        status |= gr_mat_det_cofactor(cp, mat, ctx);
        status |= gr_add(POL(1), MAT(0, 0), MAT(1, 1), ctx);
        status |= gr_neg(POL(1), POL(1), ctx);
        status |= gr_one(POL(2), ctx);
    }
    else
    {
        slong i, k, t;
        gr_ptr a, A, s;

        GR_TMP_INIT_VEC(a, n * n, ctx);
        A = GR_ENTRY(a, (n - 1) * n, sz);

        status |= _gr_vec_zero(POL(0), n + 1, ctx);
        status |= gr_neg(POL(0), MAT(0, 0), ctx);

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                status |= gr_set(GR_ENTRY(a, 0 * n + i, sz), MAT(i, t), ctx);
            }

            status |= gr_set(A, MAT(t, t), ctx);

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    s = GR_ENTRY(a, k * n + i, sz);
                    status |= _gr_vec_dot(s, NULL, 0, MAT(i, 0), GR_ENTRY(a, (k - 1) * n, sz), t + 1, ctx);
                }

                status |= gr_set(GR_ENTRY(A, k, sz), GR_ENTRY(a, k * n + t, sz), ctx);
            }

            status |= _gr_vec_dot(GR_ENTRY(A, t, sz), NULL, 0, MAT(t, 0), GR_ENTRY(a, (t - 1) * n, sz), t + 1, ctx);

            for (k = 0; k <= t; k++)
            {
                status |= _gr_vec_dot_rev(POL(k), POL(k), 1, A, cp, k, ctx);
                status |= gr_sub(POL(k), POL(k), GR_ENTRY(A, k, sz), ctx);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
            gr_swap(POL(i), POL(i - 1), ctx);

        status |= gr_one(POL(0), ctx);
        status |= _gr_poly_reverse(POL(0), POL(0), n + 1, n + 1, ctx);

        GR_TMP_CLEAR_VEC(a, n * n, ctx);
    }

#undef MAT
#undef POL

    return status;
}

int gr_mat_charpoly_berkowitz(gr_poly_t cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    int status;

    if (mat->r != mat->c)
        return GR_DOMAIN;

    gr_poly_fit_length(cp, mat->r + 1, ctx);
    _gr_poly_set_length(cp, mat->r + 1, ctx);
    status = _gr_mat_charpoly_berkowitz(cp->coeffs, mat, ctx);
    _gr_poly_normalise(cp, ctx);   /* only needed for the zero ring */
    return status;
}
