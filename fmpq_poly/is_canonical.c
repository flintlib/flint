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

int _fmpq_poly_is_canonical(const fmpz * poly, const fmpz_t den, slong len)
{
    if (len)
    {
        int ans;
        fmpz_t c;

        if (fmpz_is_zero(poly + len - 1))
            return 0;

        if (fmpz_sgn(den) < 0)
            return 0;

        fmpz_init(c);
        _fmpz_poly_content(c, poly, len);
        fmpz_gcd(c, c, den);
        ans = (*c == WORD(1));
        fmpz_clear(c);

        return ans;
    }
    else
    {
        return (*den == WORD(1));
    }
}

int fmpq_poly_is_canonical(const fmpq_poly_t poly)
{
    return _fmpq_poly_is_canonical(poly->coeffs, poly->den, poly->length);
}

