/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

void fmpz_mat_fmpz_vec_mul(
    fmpz * c,
    const fmpz * a, slong alen,
    const fmpz_mat_t B)
{
    slong len = FLINT_MIN(B->r, alen);
    fmpz_mat_struct A1, B1, C1;

    if (len <= 0)
    {
        _fmpz_vec_zero(c, B->c);
        return;
    }

    A1.entries = (fmpz *) a;
    A1.r = 1;
    A1.c = len;
    A1.stride = len;

    B1 = *B;
    B1.r = len;

    C1.entries = c;
    C1.r = 1;
    C1.c = B->c;
    C1.stride = B->c;

    fmpz_mat_mul(&C1, &A1, &B1);
}

void fmpz_mat_fmpz_vec_mul_ptr(
    fmpz * const * c,
    const fmpz * const * a, slong alen,
    const fmpz_mat_t B)
{
    slong i;
    slong len = FLINT_MIN(B->r, alen);
    slong ncols = B->c;
    fmpz * aa, * cc;
    TMP_INIT;

    len = FLINT_MAX(len, 0);

    TMP_START;

    aa = TMP_ARRAY_ALLOC(len + ncols, fmpz);
    cc = aa + len;

    /* shallow copies */
    for (i = 0; i < len; i++)
        aa[i] = *a[i];

    for (i = 0; i < ncols; i++)
        fmpz_init(cc + i);

    fmpz_mat_fmpz_vec_mul(cc, aa, len, B);

    for (i = 0; i < ncols; i++)
    {
        fmpz_swap(c[i], cc + i);
        fmpz_clear(cc + i);
    }

    TMP_END;
}
