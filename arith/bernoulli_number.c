/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void _arith_bernoulli_number(fmpz_t num, fmpz_t den, ulong n)
{
    _arith_bernoulli_number_zeta(num, den, n);
}

void arith_bernoulli_number(fmpq_t x, ulong n)
{
    _arith_bernoulli_number(fmpq_numref(x), fmpq_denref(x), n);
}
