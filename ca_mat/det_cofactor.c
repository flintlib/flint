/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

static void
_ca_mat_det_cofactor_3x3(ca_t t, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_t a, u;
    ca_init(a, ctx);
    ca_init(u, ctx);

    ca_mul   (a, ca_mat_entry(A, 1, 0), ca_mat_entry(A, 2, 1), ctx);
    ca_mul   (u, ca_mat_entry(A, 1, 1), ca_mat_entry(A, 2, 0), ctx);
    ca_sub   (a, a, u, ctx);
    ca_mul   (t, a, ca_mat_entry(A, 0, 2), ctx);

    ca_mul   (a, ca_mat_entry(A, 1, 2), ca_mat_entry(A, 2, 0), ctx);
    ca_mul   (u, ca_mat_entry(A, 1, 0), ca_mat_entry(A, 2, 2), ctx);
    ca_sub   (a, a, u, ctx);
    ca_mul   (u, a, ca_mat_entry(A, 0, 1), ctx);
    ca_add   (t, t, u, ctx);

    ca_mul   (a, ca_mat_entry(A, 1, 1), ca_mat_entry(A, 2, 2), ctx);
    ca_mul   (u, ca_mat_entry(A, 1, 2), ca_mat_entry(A, 2, 1), ctx);
    ca_sub   (a, a, u, ctx);
    ca_mul   (u, a, ca_mat_entry(A, 0, 0), ctx);
    ca_add   (t, t, u, ctx);

    ca_clear(a, ctx);
    ca_clear(u, ctx);
}

static void
ca_fmms(ca_t x, ca_t tmp, const ca_t a, const ca_t b, const ca_t c, const ca_t d, ca_ctx_t ctx)
{
    ca_mul(tmp, a, b, ctx);
    ca_mul(x, c, d, ctx);
    ca_sub(x, tmp, x, ctx);
}

static void
_ca_addmul(ca_t x, ca_t tmp, const ca_t a, const ca_t b, ca_ctx_t ctx)
{
    ca_mul(tmp, a, b, ctx);
    ca_add(x, x, tmp, ctx);
}


#define E(i,j) ca_mat_entry(A, i, j)

static void
_ca_mat_det_cofactor_4x4(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_t a, b, t;
    ca_init(a, ctx);
    ca_init(b, ctx);
    ca_init(t, ctx);

    ca_fmms(a, t, E(0,3), E(1,2), E(0,2), E(1,3), ctx);
    ca_fmms(b, t, E(2,1), E(3,0), E(2,0), E(3,1), ctx);
    ca_mul(det, a, b, ctx);

    ca_fmms(a, t, E(0,1), E(1,3), E(0,3), E(1,1), ctx);
    ca_fmms(b, t, E(2,2), E(3,0), E(2,0), E(3,2), ctx);
    _ca_addmul(det, t, a, b, ctx);

    ca_fmms(a, t, E(0,2), E(1,1), E(0,1), E(1,2), ctx);
    ca_fmms(b, t, E(2,3), E(3,0), E(2,0), E(3,3), ctx);
    _ca_addmul(det, t, a, b, ctx);

    ca_fmms(a, t, E(0,3), E(1,0), E(0,0), E(1,3), ctx);
    ca_fmms(b, t, E(2,2), E(3,1), E(2,1), E(3,2), ctx);
    _ca_addmul(det, t, a, b, ctx);

    ca_fmms(a, t, E(0,0), E(1,2), E(0,2), E(1,0), ctx);
    ca_fmms(b, t, E(2,3), E(3,1), E(2,1), E(3,3), ctx);
    _ca_addmul(det, t, a, b, ctx);

    ca_fmms(a, t, E(0,1), E(1,0), E(0,0), E(1,1), ctx);
    ca_fmms(b, t, E(2,3), E(3,2), E(2,2), E(3,3), ctx);
    _ca_addmul(det, t, a, b, ctx);

    ca_clear(a, ctx);
    ca_clear(b, ctx);
    ca_clear(t, ctx);
}

#undef E

void
ca_mat_det_cofactor(ca_t res, const ca_mat_t A, ca_ctx_t ctx)
{
    slong n;

    n = ca_mat_nrows(A);

    if (n == 0)
    {
        ca_one(res, ctx);
    }
    else if (n == 1)
    {
        ca_set(res, ca_mat_entry(A, 0, 0), ctx);
    }
    else if (n == 2)
    {
        ca_t t;
        ca_init(t, ctx);
        ca_mul(t, ca_mat_entry(A, 0, 0), ca_mat_entry(A, 1, 1), ctx);
        ca_mul(res, ca_mat_entry(A, 0, 1), ca_mat_entry(A, 1, 0), ctx);
        ca_sub(res, t, res, ctx);
        ca_clear(t, ctx);
    }
    else if (n == 3)
    {
        _ca_mat_det_cofactor_3x3(res, A, ctx);
    }
    else if (n == 4)
    {
        _ca_mat_det_cofactor_4x4(res, A, ctx);
    }
    else
    {
        flint_abort();
    }
}
