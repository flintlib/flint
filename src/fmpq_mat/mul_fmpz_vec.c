/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"


void fmpq_mat_mul_fmpz_vec(
    fmpq* c,
    const fmpq_mat_t A,
    const fmpz * b, slong blen)
{
    fmpq_t t;
    slong i, j;
    slong len = FLINT_MIN(A->c, blen);

    if (len < 1)
    {
        for (i = 0; i < A->r; i++)
            fmpq_zero(c + i);
        return;
    }

    fmpq_init(t);

    for (i = 0; i < A->r; i++)
    {
        const fmpq* Ai = A->rows[i];

        fmpq_mul_fmpz(c + i, Ai + 0, b + 0);
        for (j = 1; j < len; j++)
        {
            fmpq_mul_fmpz(t, Ai + j, b + j);
            fmpq_add(c + i, c + i, t);
        }
    }

    fmpq_clear(t);
}


void fmpq_mat_mul_fmpz_vec_ptr(
    fmpq * const * c,
    const fmpq_mat_t A,
    const fmpz * const * b, slong blen)
{
    fmpq_t t;
    slong i, j;
    slong len = FLINT_MIN(A->c, blen);

    if (len < 1)
    {
        for (i = 0; i < A->r; i++)
            fmpq_zero(c[i]);

        return;
    }

    fmpq_init(t);

    for (i = 0; i < A->r; i++)
    {
        const fmpq* Ai = A->rows[i];

        fmpq_mul_fmpz(c[i], Ai + 0, b[0]);
        for (j = 1; j < len; j++)
        {
            fmpq_mul_fmpz(t, Ai + j, b[j]);
            fmpq_add(c[i], c[i], t);
        }
    }

    fmpq_clear(t);
}

