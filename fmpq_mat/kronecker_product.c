/*
    Copyright (C) 2018 Martin Raum

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

void fmpq_mat_kronecker_product(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)
{
    slong i, j, k, l;
    slong ir, jc;

    fmpq * Aentry;

    for (i = 0, ir = 0; i < A->r; i++, ir += B->r)
    {
        for (j = 0, jc = 0; j < A->c; j++, jc += B->c)
        {
            Aentry = fmpq_mat_entry(A, i, j);
            for (k = 0; k < B->r; k++)
            {
                for (l = 0; l < B->c; l++)
                {
                    fmpq_mul(fmpq_mat_entry(C, ir + k, jc + l), Aentry, fmpq_mat_entry(B, k, l));
                }
            }
        }
    }
}
