/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_canonicalise(fmpz * poly, fmpz_t den, slong len)
{
    if (*den == WORD(1))
        return;
    
    if (*den == WORD(-1))
    {
        _fmpz_vec_neg(poly, poly, len);
        fmpz_one(den);
    }
    else if (len == 0)
    {
        fmpz_one(den);
    }
    else
    {
        fmpz_t gcd;
        fmpz_init(gcd);
        _fmpz_vec_content_chained(gcd, poly, len, den);
        if (fmpz_sgn(den) < 0)
            fmpz_neg(gcd, gcd);
        if (!fmpz_is_one(gcd))
        {
            _fmpz_vec_scalar_divexact_fmpz(poly, poly, len, gcd);
            fmpz_divexact(den, den, gcd);
        }
        fmpz_clear(gcd);
    }
}

void fmpq_poly_canonicalise(fmpq_poly_t poly)
{
    _fmpq_poly_normalise(poly);
    _fmpq_poly_canonicalise(poly->coeffs, poly->den, poly->length);
}

