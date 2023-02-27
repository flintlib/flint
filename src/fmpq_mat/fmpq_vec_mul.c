/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_vec.h"
#include "fmpq_mat.h"


void fmpq_mat_fmpq_vec_mul(
    fmpq* c,
    const fmpq* a, slong alen,
    const fmpq_mat_t B)
{
    fmpz den[1], * num;
    slong i, len = FLINT_MIN(B->r, alen);
    TMP_INIT;

    if (len < 1)
    {
        for (i = 0; i < B->c; i++)
            fmpq_zero(c + i);
        return;
    }

    TMP_START;

    fmpz_init(den);
    num = TMP_ARRAY_ALLOC(len, fmpz);
    for (i = 0; i < len; i++)
        fmpz_init(num + i);

    _fmpq_vec_get_fmpz_vec_fmpz(num, den, a, len);
    fmpq_mat_fmpz_vec_mul(c, num, len, B);

    for (i = 0; i < B->c; i++)
        fmpq_div_fmpz(c + i, c + i, den);

    fmpz_clear(den);
    for (i = 0; i < len; i++)
        fmpz_clear(num + i);

    TMP_END;
}


void fmpq_mat_fmpq_vec_mul_ptr(
    fmpq * const * c,
    const fmpq * const * a, slong alen,
    const fmpq_mat_t B)
{
    fmpz den[1], * num, ** num_ptrs;
    fmpq * ta;
    slong i, len = FLINT_MIN(B->r, alen);
    TMP_INIT;

    if (len < 1)
    {
        for (i = 0; i < B->c; i++)
            fmpq_zero(c[i]);
        return;
    }

    TMP_START;

    fmpz_init(den);
    num = TMP_ARRAY_ALLOC(len, fmpz);
    num_ptrs = TMP_ARRAY_ALLOC(len, fmpz*);
    ta = TMP_ARRAY_ALLOC(len, fmpq);

    for (i = 0; i < len; i++)
    {
        fmpz_init(num + i);
        num_ptrs[i] = num + i;
        ta[i] = *a[i];
    }

    _fmpq_vec_get_fmpz_vec_fmpz(num, den, ta, len);
    fmpq_mat_fmpz_vec_mul_ptr(c, (const fmpz * const *)num_ptrs, len, B);

    for (i = 0; i < B->c; i++)
        fmpq_div_fmpz(c[i], c[i], den);

    fmpz_clear(den);
    for (i = 0; i < len; i++)
        fmpz_clear(num + i);

    TMP_END;
}

