/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"
#include "fmpq.h"

void _arith_harmonic_number(fmpz_t num, fmpz_t den, slong n)
{
    if (n <= 0)
    {
        fmpz_zero(num);
        fmpz_one(den);
    }
    else
    {
        _fmpq_harmonic_ui(num, den, n);
    }
}

void arith_harmonic_number(fmpq_t x, slong n)
{
    _arith_harmonic_number(fmpq_numref(x), fmpq_denref(x), n);
}
