/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_poly.h"

void fmpq_poly_get_coeff_fmpq(fmpq_t x, const fmpq_poly_t poly, slong n)
{
    if (n >= poly->length)  /* Coefficient is beyond the end of poly */
    {
        fmpq_zero(x);
        return;
    }

    fmpz_set(fmpq_numref(x), poly->coeffs + n);
    fmpz_set(fmpq_denref(x), poly->den);
    fmpq_canonicalise(x);
}

void fmpq_poly_get_coeff_fmpz(fmpz_t x, const fmpq_poly_t poly, slong n)
{
    if (n >= poly->length)  /* Coefficient is beyond the end of poly */
        fmpz_zero(x);
    else
        fmpz_set(x, poly->coeffs + n);
}
