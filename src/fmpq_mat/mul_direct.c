/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"

void fmpq_mat_mul_direct(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)
{
    slong i, j, k;

    if (A == C || B == C)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_mat_mul_direct). Aliasing not implemented.\n");
    }

    if (A->c == 0)
    {
        fmpq_mat_zero(C);
        return;
    }

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            fmpq_mul(fmpq_mat_entry(C, i, j),
                     fmpq_mat_entry(A, i, 0),
                     fmpq_mat_entry(B, 0, j));

            for (k = 1; k < A->c; k++)
            {
                fmpq_addmul(fmpq_mat_entry(C, i, j),
                            fmpq_mat_entry(A, i, k),
                            fmpq_mat_entry(B, k, j));
            }
        }
    }
}
