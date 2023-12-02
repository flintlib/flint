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
fmpz_mat_randintrel(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)
{
    const slong c = mat->c, r = mat->r;

    slong i, j;

    if (c != r + 1)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_randintrel).  c != r + 1.\n");
    }

    for (i = 0; i < r; i++)
    {
        fmpz_randbits(mat->rows[i], state, bits);
        for (j = 1; j <= i; j++)
            fmpz_zero(mat->rows[i] + j);
        fmpz_one(mat->rows[i] + i + 1);
        for (j = i + 2; j < c; j++)
            fmpz_zero(mat->rows[i] + j);
    }
}
