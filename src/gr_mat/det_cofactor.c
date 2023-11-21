/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

#define E(i,j) GR_MAT_ENTRY(A, i, j, sz)

/* todo: use provided addmul/fmms */

static int
gr_fmms(gr_ptr x, gr_ptr tmp, const gr_ptr a, const gr_ptr b, const gr_ptr c, const gr_ptr d, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_mul(tmp, a, b, ctx);
    status |= gr_mul(x, c, d, ctx);
    status |= gr_sub(x, tmp, x, ctx);
    return status;
}

static int
_gr_addmul(gr_ptr x, gr_ptr tmp, const gr_ptr a, const gr_ptr b, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_mul(tmp, a, b, ctx);
    status |= gr_add(x, x, tmp, ctx);
    return status;
}

static int
_gr_mat_det_2x2(gr_ptr det, const gr_mat_t A, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr t;
    slong sz = ctx->sizeof_elem;

    GR_TMP_INIT(t, ctx);

    status |= gr_fmms(det, t, E(0,0), E(1,1), E(0,1), E(1,0), ctx);

    GR_TMP_CLEAR(t, ctx);

    return status;
}

static int
_gr_mat_det_cofactor_3x3(gr_ptr det, const gr_mat_t A, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr t, u;
    slong sz = ctx->sizeof_elem;

    GR_TMP_INIT2(t, u, ctx);

    status |= gr_fmms(t, u, E(1,0), E(2,1), E(1,1), E(2,0), ctx);
    status |= gr_mul(det, t, E(0,2), ctx);

    status |= gr_fmms(t, u, E(1,2), E(2,0), E(1,0), E(2,2), ctx);
    status |= _gr_addmul(det, u, t, E(0,1), ctx);

    status |= gr_fmms(t, u, E(1,1), E(2,2), E(1,2), E(2,1), ctx);
    status |= _gr_addmul(det, u, t, E(0,0), ctx);

    GR_TMP_CLEAR2(t, u, ctx);

    return status;
}

static int
_gr_mat_det_cofactor_4x4(gr_ptr det, const gr_mat_t A, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr a, b, t;
    slong sz = ctx->sizeof_elem;

    GR_TMP_INIT3(a, b, t, ctx);

    status |= gr_fmms(a, t, E(0,3), E(1,2), E(0,2), E(1,3), ctx);
    status |= gr_fmms(b, t, E(2,1), E(3,0), E(2,0), E(3,1), ctx);
    status |= gr_mul(det, a, b, ctx);

    status |= gr_fmms(a, t, E(0,1), E(1,3), E(0,3), E(1,1), ctx);
    status |= gr_fmms(b, t, E(2,2), E(3,0), E(2,0), E(3,2), ctx);
    status |= _gr_addmul(det, t, a, b, ctx);

    status |= gr_fmms(a, t, E(0,2), E(1,1), E(0,1), E(1,2), ctx);
    status |= gr_fmms(b, t, E(2,3), E(3,0), E(2,0), E(3,3), ctx);
    status |= _gr_addmul(det, t, a, b, ctx);

    status |= gr_fmms(a, t, E(0,3), E(1,0), E(0,0), E(1,3), ctx);
    status |= gr_fmms(b, t, E(2,2), E(3,1), E(2,1), E(3,2), ctx);
    status |= _gr_addmul(det, t, a, b, ctx);

    status |= gr_fmms(a, t, E(0,0), E(1,2), E(0,2), E(1,0), ctx);
    status |= gr_fmms(b, t, E(2,3), E(3,1), E(2,1), E(3,3), ctx);
    status |= _gr_addmul(det, t, a, b, ctx);

    status |= gr_fmms(a, t, E(0,1), E(1,0), E(0,0), E(1,1), ctx);
    status |= gr_fmms(b, t, E(2,3), E(3,2), E(2,2), E(3,3), ctx);
    status |= _gr_addmul(det, t, a, b, ctx);

    GR_TMP_CLEAR3(a, b, t, ctx);

    return status;
}

int
gr_mat_det_cofactor(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong n;

    n = A->r;

    if (n != A->c)
        return GR_DOMAIN;

    if (n == 0)
        return gr_one(res, ctx);
    else if (n == 1)
        return gr_set(res, A->rows[0], ctx);
    else if (n == 2)
        return _gr_mat_det_2x2(res, A, ctx);
    else if (n == 3)
        return _gr_mat_det_cofactor_3x3(res, A, ctx);
    else if (n == 4)
        return _gr_mat_det_cofactor_4x4(res, A, ctx);
    else
        return GR_UNABLE;
}
