/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void fmpq_div_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)
{
    fmpz_t y;

    if (fmpz_is_zero(x))
    {
        flint_printf("Exception (fmpq_div_fmpz). Division by zero.\n");
        flint_abort();
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

