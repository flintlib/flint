/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void fmpq_poly_neg(fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    if (poly1 == poly2)
    {
        _fmpz_vec_neg(poly1->coeffs, poly2->coeffs, poly2->length);
    }
    else
    {
        fmpq_poly_fit_length(poly1, poly2->length);
        _fmpz_vec_neg(poly1->coeffs, poly2->coeffs, poly2->length);
        _fmpq_poly_set_length(poly1, poly2->length);
        fmpz_set(poly1->den, poly2->den);
    }
}

