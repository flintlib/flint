/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void fmpq_div_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)
{
    fmpz_t y;

    fmpz_init(y);
    fmpz_one(y);

    _fmpq_mul(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op), fmpq_denref(op), y, x);

    fmpz_clear(y);

    if (fmpz_sgn(fmpq_denref(res)) < 0)
    {
        fmpz_neg(fmpq_numref(res), fmpq_numref(res));
        fmpz_neg(fmpq_denref(res), fmpq_denref(res));
    }
}

