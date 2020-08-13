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

int _fmpq_poly_is_monic(const fmpz * poly, const fmpz_t den, slong len)
{
    return (len > 0) && fmpz_equal(poly + (len - 1), den);
}

int fmpq_poly_is_monic(const fmpq_poly_t poly)
{
    return _fmpq_poly_is_monic(poly->coeffs, poly->den, poly->length);
}

