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
fmpz_poly_mat_mullow(fmpz_poly_mat_t C, const fmpz_poly_mat_t A,
    const fmpz_poly_mat_t B, slong len)
{
    slong ar, bc, br;
    slong i, j, k;
    fmpz_poly_t t;

    ar = A->r;
    br = B->r;
    bc = B->c;

    if (br == 0 || len < 1)
    {
        fmpz_poly_mat_zero(C);
        return;
    }

    if (C == A || C == B)
    {
        fmpz_poly_mat_t T;
        fmpz_poly_mat_init(T, ar, bc);
        fmpz_poly_mat_mullow(T, A, B, len);
        fmpz_poly_mat_swap_entrywise(C, T);
        fmpz_poly_mat_clear(T);
        return;
    }

    fmpz_poly_init(t);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            fmpz_poly_mullow(fmpz_poly_mat_entry(C, i, j),
                             fmpz_poly_mat_entry(A, i, 0),
                             fmpz_poly_mat_entry(B, 0, j), len);

            for (k = 1; k < br; k++)
            {
                fmpz_poly_mullow(t, fmpz_poly_mat_entry(A, i, k),
                                    fmpz_poly_mat_entry(B, k, j), len);

                fmpz_poly_add(fmpz_poly_mat_entry(C, i, j),
                              fmpz_poly_mat_entry(C, i, j), t);
            }
        }
    }

    fmpz_poly_clear(t);
}
