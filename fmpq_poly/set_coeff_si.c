/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"

void fmpq_poly_set_coeff_si(fmpq_poly_t poly, slong n, slong x)
{
    slong len = poly->length;
    const int replace = (n < len && !fmpz_is_zero(poly->coeffs + n));
    
    if (!replace && (x == WORD(0)))
        return;
    
    if (n + 1 > len)
    {
        fmpq_poly_fit_length(poly, n + 1);
        _fmpq_poly_set_length(poly, n + 1);
        flint_mpn_zero((mp_ptr) poly->coeffs + len, (n + 1) - len);
    }
    
    if (*poly->den == WORD(1))
    {
        fmpz_set_si(poly->coeffs + n, x);
        if (replace)
            _fmpq_poly_normalise(poly);
    }
    else
    {
        fmpz_mul_si(poly->coeffs + n, poly->den, x);
        if (replace)
            fmpq_poly_canonicalise(poly);
    }
}

