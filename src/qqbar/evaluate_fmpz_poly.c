/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "qqbar.h"

void
_qqbar_evaluate_fmpz_poly(qqbar_t res, const fmpz * poly, slong len, const qqbar_t x)
{
    fmpz_t den;
    fmpz_init(den);
    fmpz_one(den);
    _qqbar_evaluate_fmpq_poly(res, poly, den, len, x);
    fmpz_clear(den);
}

void
qqbar_evaluate_fmpz_poly(qqbar_t res, const fmpz_poly_t poly, const qqbar_t x)
{
    _qqbar_evaluate_fmpz_poly(res, poly->coeffs, poly->length, x);
}

