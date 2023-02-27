/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void _fmpq_set_ui(fmpz_t rnum, fmpz_t rden, ulong p, ulong q)
{
    if (q == 1 || p == 0)
    {
        fmpz_set_ui(rnum, p);
        fmpz_one(rden);
    }
    else
    {
        ulong r = n_gcd(p, q);

        fmpz_set_ui(rnum, p / r);
        fmpz_set_ui(rden, q / r);
    }
}

void fmpq_set_ui(fmpq_t res, ulong p, ulong q)
{
    _fmpq_set_ui(fmpq_numref(res), fmpq_denref(res), p, q);
}
