/*
    Copyright (C) 2012, 2024 Fredrik Johansson
    Copyright (C) 2015 Anubhav Srivastava

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

void
fmpz_mat_sqr(fmpz_mat_t B, const fmpz_mat_t A)
{
    slong n = A->r, ab;

    if (B == A)
    {
        fmpz_mat_t t;
        fmpz_mat_init(t, n, n);
        fmpz_mat_sqr(t, A);
        fmpz_mat_swap_entrywise(B, t);
        fmpz_mat_clear(t);
        return;
    }

    if (n <= 2)
    {
        fmpz_mat_sqr_bodrato(B, A);
        return;
    }

    ab = fmpz_mat_max_bits(A);

    if ((n == 3 && ab >= 128) || (n == 4 && ab >= 1024))
        fmpz_mat_sqr_bodrato(B, A);
    else
        fmpz_mat_mul(B, A, A);
}
