/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_approx_mul_classical(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong ar, br, bc, i, j, k;

    ar = arb_mat_nrows(A);
    br = arb_mat_nrows(B);
    bc = arb_mat_ncols(B);

    if (br == 0)
    {
        arb_mat_zero(C);
        return;
    }

    if (A == C || B == C)
    {
        arb_mat_t T;
        arb_mat_init(T, ar, bc);
        arb_mat_approx_mul_classical(T, A, B, prec);
        arb_mat_swap_entrywise(T, C);
        arb_mat_clear(T);
        return;
    }

    if (br <= 2)
    {
        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                arf_mul(arb_midref(arb_mat_entry(C, i, j)),
                          arb_midref(arb_mat_entry(A, i, 0)),
                          arb_midref(arb_mat_entry(B, 0, j)), prec, ARB_RND);

                for (k = 1; k < br; k++)
                {
                    arf_addmul(arb_midref(arb_mat_entry(C, i, j)),
                                 arb_midref(arb_mat_entry(A, i, k)),
                                 arb_midref(arb_mat_entry(B, k, j)), prec, ARB_RND);
                }
            }
        }
    }
    else
    {
        arb_ptr tmp;
        TMP_INIT;

        TMP_START;
        tmp = TMP_ALLOC(sizeof(arb_struct) * br * bc);

        for (i = 0; i < br; i++)
            for (j = 0; j < bc; j++)
                tmp[j * br + i] = *arb_mat_entry(B, i, j);

        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                arb_approx_dot(arb_mat_entry(C, i, j), NULL, 0,
                    A->rows[i], 1, tmp + j * br, 1, br, prec);
            }
        }

        TMP_END;
    }
}

void
arb_mat_approx_mul(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong cutoff;

    /* todo: detect small-integer matrices */
    if (prec <= 2 * FLINT_BITS)
        cutoff = 120;
    else if (prec <= 16 * FLINT_BITS)
        cutoff = 60;
    else
        cutoff = 40;

    if (arb_mat_nrows(A) <= cutoff || arb_mat_ncols(A) <= cutoff ||
        arb_mat_ncols(B) <= cutoff)
    {
        arb_mat_approx_mul_classical(C, A, B, prec);
    }
    else
    {
        if (arb_mat_is_exact(A) && arb_mat_is_exact(B))
        {
            arb_mat_mul(C, A, B, prec);
        }
        else
        {
            arb_mat_t AM, BM;

            if (arb_mat_is_exact(A))
            {
                arb_mat_init(BM, arb_mat_nrows(B), arb_mat_ncols(B));
                arb_mat_get_mid(BM, B);
                arb_mat_mul(C, A, BM, prec);
                arb_mat_clear(BM);
            }
            else if (arb_mat_is_exact(B))
            {
                arb_mat_init(AM, arb_mat_nrows(A), arb_mat_ncols(A));
                arb_mat_get_mid(AM, A);
                arb_mat_mul(C, AM, B, prec);
                arb_mat_clear(AM);
            }
            else
            {
                arb_mat_init(BM, arb_mat_nrows(B), arb_mat_ncols(B));
                arb_mat_get_mid(BM, B);
                arb_mat_init(AM, arb_mat_nrows(A), arb_mat_ncols(A));
                arb_mat_get_mid(AM, A);
                arb_mat_mul(C, AM, BM, prec);
                arb_mat_clear(AM);
                arb_mat_clear(BM);
            }
        }

        arb_mat_get_mid(C, C);
    }
}
