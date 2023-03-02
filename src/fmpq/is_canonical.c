/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

int _fmpq_is_canonical(const fmpz_t num, const fmpz_t den)
{
    fmpz_t u;
    int result;

    if (fmpz_is_one(den)) return 1;
    if (fmpz_sgn(den) <= 0) return 0;
    if (fmpz_is_zero(num)) return fmpz_is_one(den);

    fmpz_init(u);
    fmpz_gcd(u, num, den);
    result = fmpz_is_one(u);
    fmpz_clear(u);
    return result;
}

int fmpq_is_canonical(const fmpq_t x)
{
    return _fmpq_is_canonical(fmpq_numref(x), fmpq_denref(x));
}
