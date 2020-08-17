/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

/* sets a to gcd(a,b) and b to lcm(a,b) using temporary fmpz_t t */
static void _gcdlcm(fmpz_t t, fmpz_t a, fmpz_t b)
{
    if (fmpz_equal(a, b)) return;
    fmpz_gcd(t, a, b);
    fmpz_divexact(b, b, t);
    fmpz_mul(b, b, a);
    fmpz_set(a, t);
}

void fmpz_mat_snf_diagonal(fmpz_mat_t S, const fmpz_mat_t A)
{
    fmpz_t t;
    slong i, j, n = FLINT_MIN(A->r, A->c);

    fmpz_init(t);
    fmpz_mat_set(S, A);
    for (i = 0; i < n; i++)
        fmpz_abs(fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i, i));
    for (j = n - 1; j >= 0; j--)
    {
        for (i = 0; i < j; i++)
        {
            _gcdlcm(t, fmpz_mat_entry(S, i, i),
                    fmpz_mat_entry(S, i + 1, i + 1));
        }
    }
    fmpz_clear(t);
}
