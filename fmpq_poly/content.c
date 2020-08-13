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

void _fmpq_poly_content(fmpq_t res, const fmpz * poly, 
                        const fmpz_t den, slong len)
{
    _fmpz_poly_content(fmpq_numref(res), poly, len);
    fmpz_set(fmpq_denref(res), den);
}

void fmpq_poly_content(fmpq_t res, const fmpq_poly_t poly)
{
    _fmpq_poly_content(res, poly->coeffs, poly->den, poly->length);
}

