/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_pow_trunc(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, ulong exp,
                            slong len)
{
    slong d = fmpz_poly_mat_nrows(A);

    if (len < 1)
    {
        fmpz_poly_mat_zero(B);
    }
    else if (exp == 0 || d == 0)
    {
        fmpz_poly_mat_one(B);
    }
    else if (exp == 1)
    {
        fmpz_poly_mat_set(B, A);
        fmpz_poly_mat_truncate(B, len);
    }
    else if (exp == 2)
    {
        fmpz_poly_mat_sqrlow(B, A, len);
    }
    else if (d == 1)
    {
        fmpz_poly_pow_trunc(fmpz_poly_mat_entry(B, 0, 0),
                        fmpz_poly_mat_entry(A, 0, 0), exp, len);
    }
    else
    {
        fmpz_poly_mat_t T, U;
        slong i;

        fmpz_poly_mat_init_set(T, A);
        fmpz_poly_mat_truncate(T, len);
        fmpz_poly_mat_init(U, d, d);

        for (i = ((slong) FLINT_BIT_COUNT(exp)) - 2; i >= 0; i--)
        {
            fmpz_poly_mat_sqrlow(U, T, len);

            if (exp & (WORD(1) << i))
                fmpz_poly_mat_mullow(T, U, A, len);
            else
                fmpz_poly_mat_swap(T, U);
        }

        fmpz_poly_mat_swap(B, T);
        fmpz_poly_mat_clear(T);
        fmpz_poly_mat_clear(U);
    }
}
