/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_approx_mul_classical(acb_mat_t C, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong ar, br, bc, i, j;

    ar = acb_mat_nrows(A);
    br = acb_mat_nrows(B);
    bc = acb_mat_ncols(B);

    if (br == 0)
    {
        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                arf_zero(arb_midref(acb_realref(acb_mat_entry(C, i, j))));
                arf_zero(arb_midref(acb_imagref(acb_mat_entry(C, i, j))));
            }
        }
        return;
    }

    if (A == C || B == C)
    {
        acb_mat_t T;
        acb_mat_init(T, ar, bc);
        acb_mat_approx_mul_classical(T, A, B, prec);
        acb_mat_swap_entrywise(T, C);
        acb_mat_clear(T);
        return;
    }

    if (br == 1)
    {
        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                arf_complex_mul(
                    arb_midref(acb_realref(acb_mat_entry(C, i, j))),
                    arb_midref(acb_imagref(acb_mat_entry(C, i, j))),
                    arb_midref(acb_realref(acb_mat_entry(A, i, 0))),
                    arb_midref(acb_imagref(acb_mat_entry(A, i, 0))),
                    arb_midref(acb_realref(acb_mat_entry(B, 0, j))),
                    arb_midref(acb_imagref(acb_mat_entry(B, 0, j))), prec, ARB_RND);
            }
        }
    }
    else
    {
        acb_ptr tmp;
        TMP_INIT;

        TMP_START;
        tmp = TMP_ALLOC(sizeof(acb_struct) * br * bc);

        for (i = 0; i < br; i++)
            for (j = 0; j < bc; j++)
                tmp[j * br + i] = *acb_mat_entry(B, i, j);

        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                acb_approx_dot(acb_mat_entry(C, i, j), NULL, 0,
                    A->rows[i], 1, tmp + j * br, 1, br, prec);
            }
        }

        TMP_END;
    }
}

void
acb_mat_approx_mul(acb_mat_t C, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong cutoff;

    /* todo: detect small-integer matrices */
    if (prec <= 2 * FLINT_BITS)
        cutoff = 120;
    else if (prec <= 16 * FLINT_BITS)
        cutoff = 60;
    else
        cutoff = 40;

    if (acb_mat_nrows(A) <= cutoff || acb_mat_ncols(A) <= cutoff ||
        acb_mat_ncols(B) <= cutoff)
    {
        acb_mat_approx_mul_classical(C, A, B, prec);
    }
    else
    {
        if (acb_mat_is_exact(A) && acb_mat_is_exact(B))
        {
            acb_mat_mul(C, A, B, prec);
        }
        else
        {
            acb_mat_t AM, BM;

            if (acb_mat_is_exact(A))
            {
                acb_mat_init(BM, acb_mat_nrows(B), acb_mat_ncols(B));
                acb_mat_get_mid(BM, B);
                acb_mat_mul(C, A, BM, prec);
                acb_mat_clear(BM);
            }
            else if (acb_mat_is_exact(B))
            {
                acb_mat_init(AM, acb_mat_nrows(A), acb_mat_ncols(A));
                acb_mat_get_mid(AM, A);
                acb_mat_mul(C, AM, B, prec);
                acb_mat_clear(AM);
            }
            else
            {
                acb_mat_init(BM, acb_mat_nrows(B), acb_mat_ncols(B));
                acb_mat_get_mid(BM, B);
                acb_mat_init(AM, acb_mat_nrows(A), acb_mat_ncols(A));
                acb_mat_get_mid(AM, A);
                acb_mat_mul(C, AM, BM, prec);
                acb_mat_clear(AM);
                acb_mat_clear(BM);
            }
        }

        acb_mat_get_mid(C, C);
    }
}
