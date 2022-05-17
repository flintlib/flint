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


void fmpq_mat_mul_fmpq_vec(
    fmpq * c,
    const fmpq_mat_t A,
    const fmpq * b, slong blen)
{
    fmpz den[1], * num;
    slong i, len = FLINT_MIN(A->c, blen);
    TMP_INIT;

    if (A->r < 1 || len < 1)
    {
        for (i = 0; i < A->r; i++)
            fmpq_zero(c + i);
        return;
    }

    TMP_START;

    fmpz_init(den);
    num = TMP_ARRAY_ALLOC(len, fmpz);
    for (i = 0; i < len; i++)
        fmpz_init(num + i);

    _fmpq_vec_get_fmpz_vec_fmpz(num, den, b, len);
    fmpq_mat_mul_fmpz_vec(c, A, num, len);

    for (i = 0; i < A->r; i++)
        fmpq_div_fmpz(c + i, c + i, den);

    fmpz_clear(den);
    for (i = 0; i < len; i++)
        fmpz_clear(num + i);

    TMP_END;
}


void fmpq_mat_mul_fmpq_vec_ptr(
    fmpq * const * c,
    const fmpq_mat_t A,
    const fmpq * const * b, slong blen)
{
    fmpz den[1], * num, ** num_ptrs;
    fmpq * tb;
    slong i, len = FLINT_MIN(A->c, blen);
    TMP_INIT;

    if (len < 1)
    {
        for (i = 0; i < A->r; i++)
            fmpq_zero(c[i]);
        return;
    }

    TMP_START;

    fmpz_init(den);
    num = TMP_ARRAY_ALLOC(len, fmpz);
    num_ptrs = TMP_ARRAY_ALLOC(len, fmpz*);
    tb = TMP_ARRAY_ALLOC(len, fmpq);

    for (i = 0; i < len; i++)
    {
        fmpz_init(num + i);
        num_ptrs[i] = num + i;
        tb[i] = *b[i];
    }

    _fmpq_vec_get_fmpz_vec_fmpz(num, den, tb, len);
    fmpq_mat_mul_fmpz_vec_ptr(c, A, (const fmpz * const *)num_ptrs, len);

    for (i = 0; i < A->r; i++)
        fmpq_div_fmpz(c[i], c[i], den);

    fmpz_clear(den);
    for (i = 0; i < len; i++)
        fmpz_clear(num + i);

    TMP_END;
}

