/*
    Copyright (C) 2011, 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
_fmpq_mul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den,
            const fmpz_t op2num, const fmpz_t op2den)
{
    if (!COEFF_IS_MPZ(*op1num) && !COEFF_IS_MPZ(*op1den) && !COEFF_IS_MPZ(*op2num) && !COEFF_IS_MPZ(*op2den))
    {
        _fmpq_mul_small(rnum, rden, *op1num, *op1den, *op2num, *op2den);
        return;
    }

    /* Common special cases: squaring, same denominator (e.g. both integers) */
    if (((op1num == op2num) && (op1den == op2den)) ||
         fmpz_equal(op1den, op2den))
    {
        fmpz_mul(rnum, op1num, op2num);
        fmpz_mul(rden, op1den, op2den);
    }
    /* Exactly one argument is an integer */
    else if (fmpz_is_one(op1den))
    {
        fmpz_t t, x;
        fmpz_init(t);
        fmpz_init(x);

        fmpz_gcd(t, op1num, op2den);
        fmpz_divexact(x, op1num, t);
        fmpz_mul(rnum, x, op2num);
        fmpz_divexact(t, op2den, t);
        fmpz_mul(rden, op1den, t);

        fmpz_clear(t);
        fmpz_clear(x);
    }
    else if (fmpz_is_one(op2den))
    {
        fmpz_t t, x;
        fmpz_init(t);
        fmpz_init(x);

        fmpz_gcd(t, op2num, op1den);
        fmpz_divexact(x, op2num, t);
        fmpz_mul(rnum, x, op1num);
        fmpz_divexact(t, op1den, t);
        fmpz_mul(rden, op2den, t);

        fmpz_clear(t);
        fmpz_clear(x);
    }
    else
    {
        fmpz_t t, u, x, y;

        fmpz_init(t);
        fmpz_init(u);
        fmpz_init(x);
        fmpz_init(y);

        fmpz_gcd(t, op1num, op2den);
        fmpz_gcd(u, op1den, op2num);

        fmpz_divexact(x, op1num, t);
        fmpz_divexact(y, op2num, u);

        fmpz_mul(rnum, x, y);

        fmpz_divexact(x, op1den, u);
        fmpz_divexact(y, op2den, t);

        fmpz_mul(rden, x, y);

        fmpz_clear(t);
        fmpz_clear(u);
        fmpz_clear(x);
        fmpz_clear(y);
    }
}


void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
{
    _fmpq_mul(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1),
              fmpq_numref(op2), fmpq_denref(op2));
}
