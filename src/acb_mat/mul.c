/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

static slong
acb_mat_bits(const acb_mat_t A)
{
    slong b, t, i, ar, ac;

    ar = acb_mat_nrows(A);
    ac = acb_mat_ncols(A);

    b = 0;
    for (i = 0; i < ar; i++)
    {
        t = _arb_vec_bits((arb_srcptr) A->rows[i], 2 * ac);
        b = FLINT_MAX(b, t);
    }

    return b;
}

int acb_mat_is_lagom(const acb_mat_t A)
{
    slong i, j, M, N;

    M = acb_mat_nrows(A);
    N = acb_mat_ncols(A);

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (!ARB_IS_LAGOM(acb_realref(acb_mat_entry(A, i, j))) ||
                !ARB_IS_LAGOM(acb_imagref(acb_mat_entry(A, i, j))))
                return 0;
        }
    }

    return 1;
}

void
acb_mat_mul(acb_mat_t C, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong ar, ac, br, bc, n, abits, bbits, bits;

    ar = acb_mat_nrows(A);
    ac = acb_mat_ncols(A);
    br = acb_mat_nrows(B);
    bc = acb_mat_ncols(B);

    if (ac != br || ar != acb_mat_nrows(C) || bc != acb_mat_ncols(C))
    {
        flint_throw(FLINT_ERROR, "acb_mat_mul: incompatible dimensions\n");
    }

    n = FLINT_MIN(ar, ac);
    n = FLINT_MIN(ac, bc);

    if (n >= 20)
    {
        abits = acb_mat_bits(A);
        bbits = acb_mat_bits(B);

        bits = FLINT_MIN(prec, FLINT_MAX(abits, bbits));

        if (bits < 8000 && n >= 5 + bits / 64)
        {
            acb_mat_mul_reorder(C, A, B, prec);
            return;
        }
    }

    if (flint_get_num_threads() > 1 &&
        ((double) ar * (double) ac * (double) bc * (double) prec > 100000))
    {
        acb_mat_mul_threaded(C, A, B, prec);
    }
    else
    {
        acb_mat_mul_classical(C, A, B, prec);
    }
}

