/*
    Copyright (C) 2005-2009 Damien Stehle
    Copyright (C) 2007 David Cade
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

void
fmpz_mat_randsimdioph(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits, flint_bitcnt_t bits2)
{
    const slong c = mat->c, r = mat->r;

    slong i, j;

    if (c != r)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_randsimdioph). Ill-formed matrix.\n");
    }

    fmpz_one(mat->rows[0]);
    fmpz_mul_2exp(mat->rows[0], mat->rows[0], bits2);
    for (j = 1; j < c; j++)
        fmpz_randbits(mat->rows[0] + j, state, bits);
    for (i = 1; i < r; i++)
    {
        for (j = 0; j < i; j++)
            fmpz_zero(mat->rows[i] + j);
        fmpz_one(mat->rows[i] + i);
        fmpz_mul_2exp(mat->rows[i] + i, mat->rows[i] + i, bits);
        for (j = i + 1; j < c; j++)
            fmpz_zero(mat->rows[i] + j);
    }
}
