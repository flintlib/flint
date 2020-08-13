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
#include "fmpq_poly.h"

void fmpq_poly_set_mpq(fmpq_poly_t poly, const mpq_t x)
{
    fmpq_poly_fit_length(poly, 1);
    fmpz_set_mpz(poly->coeffs, mpq_numref(x));
    fmpz_set_mpz(poly->den, mpq_denref(x));
    _fmpq_poly_set_length(poly, 1);
    _fmpq_poly_normalise(poly);
}

