/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
_qqbar_get_fmpq(fmpz_t num, fmpz_t den, const qqbar_t x)
{
    if (qqbar_degree(x) != 1)
    {
        flint_throw(FLINT_ERROR, "_qqbar_get_fmpq: not a rational value\n");
    }

    fmpz_neg(num, QQBAR_COEFFS(x));
    fmpz_set(den, QQBAR_COEFFS(x) + 1);
}

void
qqbar_get_fmpq(fmpq_t res, const qqbar_t x)
{
    _qqbar_get_fmpq(fmpq_numref(res), fmpq_denref(res), x);
}
