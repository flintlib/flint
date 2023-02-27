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
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void _fmpq_poly_make_monic(fmpz * rpoly, fmpz_t rden, 
                      const fmpz * poly, const fmpz_t den, slong len)
{
    _fmpz_poly_primitive_part(rpoly, poly, len);
    fmpz_set(rden, rpoly + (len - 1));
}

void fmpq_poly_make_monic(fmpq_poly_t res, const fmpq_poly_t poly)
{
    slong len = poly->length;
    if (len == 0)
    {
        fmpq_poly_zero(res);
        return;
    }
    
    fmpq_poly_fit_length(res, len);
    _fmpq_poly_set_length(res, len);
    
    _fmpq_poly_make_monic(res->coeffs, res->den, poly->coeffs, poly->den, len);
}

