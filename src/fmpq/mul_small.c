/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
_fmpq_mul_small(fmpz_t rnum, fmpz_t rden, slong op1num, ulong op1den, slong op2num, ulong op2den)
{
    mp_limb_t hi, lo, denhi, denlo;
    int neg;

    if (op1num == 0 || op2num == 0)
    {
        fmpz_zero(rnum);
        fmpz_one(rden);
        return;
    }

    neg = 0;
    if (op1num < 0)
    {
        op1num = -op1num;
        neg = 1;
    }

    if (op2num < 0)
    {
        op2num = -op2num;
        neg = !neg;
    }

    if (op1den == op2den)
    {
        umul_ppmm(hi, lo, op1num, op2num);
        umul_ppmm(denhi, denlo, op1den, op2den);
    }
    else if (op1den == 1)
    {
        ulong t, x;
        t = n_gcd(op1num, op2den);
        x = op1num / t;
        t = op2den / t;
        umul_ppmm(hi, lo, x, op2num);
        umul_ppmm(denhi, denlo, op1den, t);
    }
    else if (op2den == 1)
    {
        ulong t, x;
        t = n_gcd(op2num, op1den);
        x = op2num / t;
        t = op1den / t;
        umul_ppmm(hi, lo, x, op1num);
        umul_ppmm(denhi, denlo, op2den, t);
    }
    else
    {
        ulong t, u, x, y;
        t = n_gcd(op1num, op2den);
        u = n_gcd(op1den, op2num);
        x = op1num / t;
        y = op2num / u;
        umul_ppmm(hi, lo, x, y);
        x = op1den / u;
        y = op2den / t;
        umul_ppmm(denhi, denlo, x, y);
    }

    if (neg)
        fmpz_neg_uiui(rnum, hi, lo);
    else
        fmpz_set_uiui(rnum, hi, lo);
    fmpz_set_uiui(rden, denhi, denlo);
}
