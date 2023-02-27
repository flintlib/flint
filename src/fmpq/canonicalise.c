/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void _fmpq_canonicalise(fmpz_t num, fmpz_t den)
{
    fmpz_t u;

    if (fmpz_is_one(den))
        return;

    if (fmpz_is_zero(num))
    {
        fmpz_one(den);
        return;
    }

    fmpz_init(u);
    fmpz_gcd(u, num, den);

    if (!fmpz_is_one(u))
    {
        fmpz_divexact(num, num, u);
        fmpz_divexact(den, den, u);
    }

    fmpz_clear(u);

    if (fmpz_sgn(den) < 0)
    {
        fmpz_neg(num, num);
        fmpz_neg(den, den);
    }
}

void fmpq_canonicalise(fmpq_t res)
{
    _fmpq_canonicalise(fmpq_numref(res), fmpq_denref(res));
}
