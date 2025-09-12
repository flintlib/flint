/*
    Copyright (C) 2025 Kacper Proniewski

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_poly.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void _fmpq_poly_discriminant(
    fmpq_t res, const fmpz* poly, const fmpz_t den, slong len)
{
    fmpz_t rnum, rden;
    fmpz_init(rnum);
    fmpz_init(rden);

    _fmpz_poly_discriminant(rnum, poly, len);

    // according to the formula disc(c*f)=c^(2n-2)disc(f) where n=len-1 is the degree
    fmpz_pow_ui(rden, den, (ulong)2 * len - 4);

    fmpq_set_fmpz_frac(res, rnum, rden);

    fmpz_clear(rnum);
    fmpz_clear(rden);
}

void fmpq_poly_discriminant(fmpq_t res, const fmpq_poly_t poly)
{
    slong len = poly->length;

    if (len <= 1)
        fmpq_zero(res);
    else
        _fmpq_poly_discriminant(res, poly->coeffs, poly->den, len);
}
