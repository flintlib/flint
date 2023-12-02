/*
    Copyright (C) 2012, 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_mul_classical(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong ar, ac, br, bc, i, j, k;

    if (A == B && (arb_mat_nrows(A) <= 2 ||
        (prec >= 1024 && arb_mat_nrows(A) < 8)))
    {
        arb_mat_sqr_classical(C, A, prec);
        return;
    }

    ar = arb_mat_nrows(A);
    ac = arb_mat_ncols(A);
    br = arb_mat_nrows(B);
    bc = arb_mat_ncols(B);

    if (ac != br || ar != arb_mat_nrows(C) || bc != arb_mat_ncols(C))
    {
        flint_throw(FLINT_ERROR, "arb_mat_mul: incompatible dimensions\n");
    }

    if (br == 0)
    {
        arb_mat_zero(C);
        return;
    }

    if (A == C || B == C)
    {
        arb_mat_t T;
        arb_mat_init(T, ar, bc);
        arb_mat_mul_classical(T, A, B, prec);
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
                /* todo: efficient fmma code */
                arb_mul(arb_mat_entry(C, i, j),
                          arb_mat_entry(A, i, 0),
                          arb_mat_entry(B, 0, j), prec);

                for (k = 1; k < br; k++)
                {
                    arb_addmul(arb_mat_entry(C, i, j),
                                 arb_mat_entry(A, i, k),
                                 arb_mat_entry(B, k, j), prec);
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
                arb_dot(arb_mat_entry(C, i, j), NULL, 0,
                    A->rows[i], 1, tmp + j * br, 1, br, prec);
            }
        }

        TMP_END;
    }
}
