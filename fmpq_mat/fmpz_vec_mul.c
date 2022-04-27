/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"


void fmpq_mat_fmpz_vec_mul(
    fmpq * c,
    const fmpz * a, slong alen,
    const fmpq_mat_t B)
{
    fmpq_t t;
    slong i, j;
    slong len = FLINT_MIN(B->r, alen);

    if (len < 1)
    {
        for (i = 0; i < B->c; i++)
            fmpq_zero(c + i);

        return;
    }

    fmpq_init(t);

    for (i = 0; i < B->c; i++)
    {
        fmpq_mul_fmpz(c + i, fmpq_mat_entry(B, 0, i), a + 0);
        for (j = 1; j < len; j++)
        {
            fmpq_mul_fmpz(t, fmpq_mat_entry(B, j, i), a + j);
            fmpq_add(c + i, c + i, t);
        }
    }

    fmpq_clear(t);
}


void fmpq_mat_fmpz_vec_mul_ptr(
    fmpq * const * c,
    const fmpz * const * a, slong alen,
    const fmpq_mat_t B)
{
    fmpq_t t;
    slong i, j;
    slong len = FLINT_MIN(B->r, alen);

    if (len < 1)
    {
        for (i = 0; i < B->c; i++)
            fmpq_zero(c[i]);

        return;
    }

    fmpq_init(t);

    for (i = 0; i < B->c; i++)
    {
        fmpq_mul_fmpz(c[i], fmpq_mat_entry(B, 0, i), a[0]);
        for (j = 1; j < len; j++)
        {
            fmpq_mul_fmpz(t, fmpq_mat_entry(B, j, i), a[j]);
            fmpq_add(c[i], c[i], t);
        }
    }

    fmpq_clear(t);
}

