/*
    Copyright (C) 2017 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_squarefree_part(fmpz_poly_t res, const fmpz_poly_t poly)
{

    if (poly->length == 0)
    {
        fmpz_poly_zero(res);
    }
    else if (poly->length == 1)
    {
        fmpz_poly_one(res);
    }
    else
    {
        fmpz_poly_t der, gcd;

        fmpz_poly_init(der);
        fmpz_poly_init(gcd);

        fmpz_poly_derivative(der, poly);
        fmpz_poly_gcd(gcd, poly, der);
        fmpz_poly_divexact(res, poly, gcd);

        if (res->length && fmpz_sgn(res->coeffs + res->length - 1) < 0)
            fmpz_poly_neg(res, res);

        fmpz_poly_clear(der);
        fmpz_poly_clear(gcd);
    }
}
