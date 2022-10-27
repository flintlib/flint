/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"


int _fmpq_cmp_fmpz(const fmpz_t p, const fmpz_t q, const fmpz_t r)
{
    fmpz one = 1;
    return _fmpq_cmp(p, q, r, &one);
}


int fmpq_cmp_fmpz(const fmpq_t x, const fmpz_t y)
{
    return _fmpq_cmp_fmpz(fmpq_numref(x), fmpq_denref(x), y);
}

