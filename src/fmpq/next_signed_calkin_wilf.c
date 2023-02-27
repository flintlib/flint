/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
_fmpq_next_signed_calkin_wilf(fmpz_t rnum, fmpz_t rden,
    const fmpz_t num, const fmpz_t den)
{
    if (fmpz_sgn(num) > 0)
    {
        fmpz_neg(rnum, num);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_neg(rnum, num);
        _fmpq_next_calkin_wilf(rnum, rden, rnum, den);
    }
}

void
fmpq_next_signed_calkin_wilf(fmpq_t res, const fmpq_t x)
{
    _fmpq_next_signed_calkin_wilf(fmpq_numref(res), fmpq_denref(res),
        fmpq_numref(x), fmpq_denref(x));
}
