/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"

void
fmpq_mat_gso(fmpq_mat_t B, const fmpq_mat_t A)
{
    slong i, j, k;
    fmpq_t num, den, mu;

    if (B->r != A->r || B->c != A->c)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_mat_gso). Incompatible dimensions.\n");
    }

    if (B == A)
    {
        fmpq_mat_t t;
        fmpq_mat_init(t, B->r, B->c);
        fmpq_mat_gso(t, A);
        fmpq_mat_swap_entrywise(B, t);
        fmpq_mat_clear(t);
        return;
    }

    if (A->r == 0)
    {
        return;
    }

    fmpq_init(num);
    fmpq_init(den);
    fmpq_init(mu);

    for (i = 0; i < A->c; i++)
    {
        for (j = 0; j < A->r; j++)
        {
            fmpq_set(fmpq_mat_entry(B, j, i), fmpq_mat_entry(A, j, i));
        }

        for (j = 0; j < i; j++)
        {
            fmpq_mul(num, fmpq_mat_entry(A, 0, i), fmpq_mat_entry(B, 0, j));

            for (k = 1; k < A->r; k++)
            {
                fmpq_addmul(num,
                            fmpq_mat_entry(A, k, i), fmpq_mat_entry(B, k, j));
            }

            fmpq_mul(den, fmpq_mat_entry(B, 0, j), fmpq_mat_entry(B, 0, j));

            for (k = 1; k < A->r; k++)
            {
                fmpq_addmul(den,
                            fmpq_mat_entry(B, k, j), fmpq_mat_entry(B, k, j));
            }

            if (!fmpq_is_zero(den))
            {
                fmpq_div(mu, num, den);

                for (k = 0; k < A->r; k++)
                {
                    fmpq_submul(fmpq_mat_entry(B, k, i),
                                mu, fmpq_mat_entry(B, k, j));
                }
            }
        }
    }

    fmpq_clear(num);
    fmpq_clear(den);
    fmpq_clear(mu);
}
