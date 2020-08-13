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

static void
binary_splitting(fmpz_poly_mat_t P, fmpz_poly_mat_t * const factors,
                                                            slong n1, slong n2)
{
    if (n2 - n1 <= 0)
    {
        fmpz_poly_mat_one(P);
    }
    else if (n2 - n1 == 1)
    {
        fmpz_poly_mat_set(P, factors[n1]);
    }
    else if (n2 - n1 == 2)
    {
        fmpz_poly_mat_mul(P, factors[n1], factors[n1 + 1]);
    }
    else
    {
        fmpz_poly_mat_t P1, P2;
        slong m = (n1 + n2) / 2;

        fmpz_poly_mat_init(P1, P->r, P->c);
        fmpz_poly_mat_init(P2, P->r, P->c);

        binary_splitting(P1, factors, n1, m);
        binary_splitting(P2, factors, m, n2);

        fmpz_poly_mat_mul(P, P1, P2);

        fmpz_poly_mat_clear(P1);
        fmpz_poly_mat_clear(P2);
    }
}

void
fmpz_poly_mat_prod(fmpz_poly_mat_t res,
                        fmpz_poly_mat_t * const factors, slong n)
{
    binary_splitting(res, factors, 0, n);
}
