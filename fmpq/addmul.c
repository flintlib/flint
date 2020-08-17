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
_fmpq_addmul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den,
            const fmpz_t op2num, const fmpz_t op2den)
{
    fmpz_t u, v;

    fmpz_init(u);
    fmpz_init(v);

    _fmpq_mul(u, v, op1num, op1den, op2num, op2den);
    _fmpq_add(rnum, rden, rnum, rden, u, v);

    fmpz_clear(u);
    fmpz_clear(v);
}

void fmpq_addmul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
{
    _fmpq_addmul(fmpq_numref(res), fmpq_denref(res),
                 fmpq_numref(op1), fmpq_denref(op1),
                 fmpq_numref(op2), fmpq_denref(op2));
}
