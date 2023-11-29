/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
_fmpq_div(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den,
                                    const fmpz_t op2num, const fmpz_t op2den)
{
    fmpz_t t, u;

    if (!COEFF_IS_MPZ(*op1num) && !COEFF_IS_MPZ(*op1den) && !COEFF_IS_MPZ(*op2num) && !COEFF_IS_MPZ(*op2den))
    {
        if (*op2num > 0)
            _fmpq_mul_small(rnum, rden, *op1num, *op1den, *op2den, *op2num);
        else
            _fmpq_mul_small(rnum, rden, *op1num, *op1den, -(*op2den), -(*op2num));
        return;
    }

    fmpz_init(t);
    fmpz_init(u);
    fmpz_set(t, op2den);
    fmpz_set(u, op2num);

    _fmpq_mul(rnum, rden, op1num, op1den, t, u);

    fmpz_clear(t);
    fmpz_clear(u);

    if (fmpz_sgn(rden) < 0)
    {
        fmpz_neg(rnum, rnum);
        fmpz_neg(rden, rden);
    }
}

void fmpq_div(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
{
    if (fmpq_is_zero(op2))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_div). Division by zero.\n");
    }

    _fmpq_div(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1),
              fmpq_numref(op2), fmpq_denref(op2));
}

void fmpq_div_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)
{
    fmpz_t y;

    if (fmpz_is_zero(x))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_div_fmpz). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(*fmpq_numref(op)) && !COEFF_IS_MPZ(*fmpq_denref(op)) && !COEFF_IS_MPZ(*x))
    {
        if (*x >= 0)
            _fmpq_mul_small(fmpq_numref(res), fmpq_denref(res), *fmpq_numref(op), *fmpq_denref(op), 1, *x);
        else
            _fmpq_mul_small(fmpq_numref(res), fmpq_denref(res), *fmpq_numref(op), *fmpq_denref(op), -WORD(1), -(*x));
        return;
    }

    *y = 1;
    _fmpq_mul(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op), fmpq_denref(op), y, x);

    if (fmpz_sgn(fmpq_denref(res)) < 0)
    {
        fmpz_neg(fmpq_numref(res), fmpq_numref(res));
        fmpz_neg(fmpq_denref(res), fmpq_denref(res));
    }
}
