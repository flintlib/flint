/*
    Copyright (C) 2012, 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_mul_classical(acb_mat_t C, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong ar, ac, br, bc, i, j, k;

    if (A == B && (acb_mat_nrows(A) <= 2 ||
        (prec >= 1024 && acb_mat_nrows(A) < 8)))
    {
        acb_mat_sqr_classical(C, A, prec);
        return;
    }

    ar = acb_mat_nrows(A);
    ac = acb_mat_ncols(A);
    br = acb_mat_nrows(B);
    bc = acb_mat_ncols(B);

    if (ac != br || ar != acb_mat_nrows(C) || bc != acb_mat_ncols(C))
    {
        flint_throw(FLINT_ERROR, "acb_mat_mul: incompatible dimensions\n");
    }

    if (br == 0)
    {
        acb_mat_zero(C);
        return;
    }

    if (A == C || B == C)
    {
        acb_mat_t T;
        acb_mat_init(T, ar, bc);
        acb_mat_mul_classical(T, A, B, prec);
        acb_mat_swap_entrywise(T, C);
        acb_mat_clear(T);
        return;
    }

    if (br <= 2)
    {
        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                /* todo: efficient fmma code */
                acb_mul(acb_mat_entry(C, i, j),
                          acb_mat_entry(A, i, 0),
                          acb_mat_entry(B, 0, j), prec);

                for (k = 1; k < br; k++)
                {
                    acb_addmul(acb_mat_entry(C, i, j),
                                 acb_mat_entry(A, i, k),
                                 acb_mat_entry(B, k, j), prec);
                }
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
                acb_dot(acb_mat_entry(C, i, j), NULL, 0,
                    A->rows[i], 1, tmp + j * br, 1, br, prec);
            }
        }

        TMP_END;
    }
}
