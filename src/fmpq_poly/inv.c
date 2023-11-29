/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq_poly.h"

void fmpq_poly_inv(fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    if (poly2->length != 1)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_inv). poly2 is not invertible.\n");
    }

    if (poly1 == poly2)
    {
        fmpz_swap(poly1->coeffs, poly1->den);
        if (fmpz_sgn(poly1->den) < 0)
        {
            fmpz_neg(poly1->coeffs, poly1->coeffs);
            fmpz_neg(poly1->den, poly1->den);
        }
    }
    else
    {
        fmpq_poly_fit_length(poly1, 1);
        if (fmpz_sgn(poly2->coeffs) > 0)
        {
            fmpz_set(poly1->coeffs, poly2->den);
            fmpz_set(poly1->den, poly2->coeffs);
        }
        else
        {
            fmpz_neg(poly1->coeffs, poly2->den);
            fmpz_neg(poly1->den, poly2->coeffs);
        }
        _fmpq_poly_set_length(poly1, 1);
    }
}

