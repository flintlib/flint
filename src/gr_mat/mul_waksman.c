/*
    Copyright (C) 2024 Éric Schost
    Copyright (C) 2024 Vincent Neiger
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"
#include "gr_vec.h"

/* todo: division by two should be divexact by two */
/* todo: avoid redundant additions 0 + ... */

int gr_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong m, n, p;
    slong sz = ctx->sizeof_elem;

    m = A->r;
    n = A->c;
    p = B->c;

    if (m == 0 || n == 0 || p == 0)
    {
        return gr_mat_zero(C, ctx);
    }

    if (n != B->r || m != C->r || p != C->c)
    {
        return GR_DOMAIN;
    }

    if (A == C || B == C)
    {
        gr_mat_t T;
        gr_mat_init(T, m, p, ctx);
        status |= gr_mat_mul_waksman(T, A, B, ctx);
        status |= gr_mat_swap_entrywise(T, C, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

    slong i, l, j, k;

    gr_ptr tmp, Crow, Ccol, val0, val1, val2, crow;

    GR_TMP_INIT_VEC(tmp, p + m + 4, ctx);

    Crow = tmp;
    Ccol = GR_ENTRY(Crow, p, sz);
    val0 = GR_ENTRY(Ccol, m, sz);
    val1 = GR_ENTRY(val0, 1, sz);
    val2 = GR_ENTRY(val1, 1, sz);
    crow = GR_ENTRY(val2, 1, sz);

    slong np = n >> 1;

    for (i = 0; i < m; i++)
        status |= _gr_vec_zero(GR_MAT_ENTRY(C, i, 0, sz), p, ctx);

    for (j = 1; j <= np; j++)
    {
        slong j2 = (j << 1) - 1;
    
        for (k = 0; k < p; k++)
        {
            status |= gr_add(val1, GR_MAT_ENTRY(A, 0, j2 - 1, sz), GR_MAT_ENTRY(B, j2, k, sz), ctx);
            status |= gr_add(val2, GR_MAT_ENTRY(A, 0, j2, sz), GR_MAT_ENTRY(B, j2 - 1, k, sz), ctx);
            status |= gr_addmul(GR_MAT_ENTRY(C, 0, k, sz), val1, val2, ctx);

            status |= gr_sub(val1, GR_MAT_ENTRY(A, 0, j2 - 1, sz), GR_MAT_ENTRY(B, j2, k, sz), ctx);
            status |= gr_sub(val2, GR_MAT_ENTRY(A, 0, j2, sz), GR_MAT_ENTRY(B, j2 - 1, k, sz), ctx);
            status |= gr_addmul(GR_ENTRY(Crow, k, sz), val1, val2, ctx);
        }

        for (l = 1; l < m; l++)
        {
            status |= gr_add(val1, GR_MAT_ENTRY(A, l, j2 - 1, sz), GR_MAT_ENTRY(B, j2, 0, sz), ctx);
            status |= gr_add(val2, GR_MAT_ENTRY(A, l, j2, sz), GR_MAT_ENTRY(B, j2 - 1, 0, sz), ctx);
            status |= gr_addmul(GR_MAT_ENTRY(C, l, 0, sz), val1, val2, ctx);
      
            status |= gr_sub(val1, GR_MAT_ENTRY(A, l, j2 - 1, sz), GR_MAT_ENTRY(B, j2, 0, sz), ctx);
            status |= gr_sub(val2, GR_MAT_ENTRY(A, l, j2, sz), GR_MAT_ENTRY(B, j2 - 1, 0, sz), ctx);
            status |= gr_addmul(GR_ENTRY(Ccol, l, sz), val1, val2, ctx);
        }

        for (k = 1; k < p; k++)
        {
            for (l = 1; l < m; l++)
            {
                status |= gr_add(val1, GR_MAT_ENTRY(A, l, j2 - 1, sz), GR_MAT_ENTRY(B, j2, k, sz), ctx);
                status |= gr_add(val2, GR_MAT_ENTRY(A, l, j2, sz), GR_MAT_ENTRY(B, j2 - 1, k, sz), ctx);
                status |= gr_addmul(GR_MAT_ENTRY(C, l, k, sz), val1, val2, ctx);
            }
        }
    }

    for (l = 1; l < m; l++)
    {
        status |= gr_add(val1, GR_ENTRY(Ccol, l, sz), GR_MAT_ENTRY(C, l, 0, sz), ctx);
        status |= gr_mul_2exp_si(GR_ENTRY(Ccol, l, sz), val1, -1, ctx); 
        status |= gr_sub(GR_MAT_ENTRY(C, l, 0, sz), GR_MAT_ENTRY(C, l, 0, sz), GR_ENTRY(Ccol, l, sz), ctx);
    }

    status |= gr_add(val1, Crow, GR_MAT_ENTRY(C, 0, 0, sz), ctx);
    status |= gr_mul_2exp_si(val0, val1, -1, ctx);
    status |= gr_sub(GR_MAT_ENTRY(C, 0, 0, sz), GR_MAT_ENTRY(C, 0, 0, sz), val0, ctx);

    for (k = 1; k < p; k++)
    {
        status |= gr_add(crow, GR_ENTRY(Crow, k, sz), GR_MAT_ENTRY(C, 0, k, sz), ctx);
        status |= gr_mul_2exp_si(val1, crow, -1, ctx); 
        status |= gr_sub(GR_MAT_ENTRY(C, 0, k, sz), GR_MAT_ENTRY(C, 0, k, sz), val1, ctx);
        status |= gr_sub(crow, val1, val0, ctx);

        for (l = 1; l < m; l++)
        {
            status |= gr_sub(val2, GR_MAT_ENTRY(C, l, k, sz), crow, ctx);
            status |= gr_sub(GR_MAT_ENTRY(C, l, k, sz), val2, GR_ENTRY(Ccol, l, sz), ctx);
        }
    }

    if ((n & 1) == 1)
        for (l = 0; l < m; l++)
            for (k = 0; k < p; k++)
                status |= gr_addmul(GR_MAT_ENTRY(C, l, k, sz),
                    GR_MAT_ENTRY(A, l, n - 1, sz), GR_MAT_ENTRY(B, n - 1, k, sz), ctx);

    GR_TMP_CLEAR_VEC(tmp, p + m + 4, ctx);

    return status;
}

