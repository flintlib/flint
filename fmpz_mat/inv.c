/*
    Copyright (C) 2010-2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

static void
_fmpz_mat_inv_2x2(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
{
    fmpz_fmms(den, fmpz_mat_entry(A, 0, 0), fmpz_mat_entry(A, 1, 1),
                   fmpz_mat_entry(A, 0, 1), fmpz_mat_entry(A, 1, 0));

    fmpz_neg(fmpz_mat_entry(B, 0, 1), fmpz_mat_entry(A, 0, 1));
    fmpz_neg(fmpz_mat_entry(B, 1, 0), fmpz_mat_entry(A, 1, 0));

    if (A == B)
    {
        fmpz_swap(fmpz_mat_entry(B, 0, 0), fmpz_mat_entry(B, 1, 1));
    }
    else
    {
        fmpz_set(fmpz_mat_entry(B, 0, 0), fmpz_mat_entry(A, 1, 1));
        fmpz_set(fmpz_mat_entry(B, 1, 1), fmpz_mat_entry(A, 0, 0));
    }
}

int
fmpz_mat_inv(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
{
    slong dim = A->r;

    if (dim == 0)
    {
        fmpz_one(den);
        return 1;
    }
    else if (dim == 1)
    {
        fmpz_set(den, fmpz_mat_entry(A, 0, 0));
        fmpz_one(fmpz_mat_entry(B, 0, 0));
        return !fmpz_is_zero(den);
    }
    else if (dim == 2)
    {
        _fmpz_mat_inv_2x2(B, den, A);
        return !fmpz_is_zero(den);
    }
    else
    {
        fmpz_mat_t I;
        slong i;
        int success;

        fmpz_mat_init(I, dim, dim);
        for (i = 0; i < dim; i++)
            fmpz_one(fmpz_mat_entry(I, i, i));
        success = fmpz_mat_solve(B, den, A, I);
        fmpz_mat_clear(I);
        return success;
    }
}
