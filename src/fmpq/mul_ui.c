/*
    Copyright (C) 2020 William Hart
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

static ulong _fmpz_gcd_ui(const fmpz_t g, ulong h)
{
    if (!COEFF_IS_MPZ(*g))
        return n_gcd(FLINT_ABS(*g), h);
    else
        return n_gcd(flint_mpz_fdiv_ui(COEFF_TO_PTR(*g), h), h);
}

void
_fmpq_mul_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            ulong r)
{
    if (r == 0 || fmpz_is_zero(p))
    {
        fmpz_zero(rnum);
        fmpz_one(rden);
    }
    else if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && r <= COEFF_MAX)
    {
        _fmpq_mul_small(rnum, rden, *p, *q, r, 1);
    }
    else if (r == 1)
    {
        fmpz_set(rnum, p);
        fmpz_set(rden, q);
    }
    else
    {
        ulong g = _fmpz_gcd_ui(q, r);

        if (g == 1)
        {
            fmpz_set(rden, q);
            fmpz_mul_ui(rnum, p, r);
        }
        else
        {
            fmpz_mul_ui(rnum, p, r / g);
            fmpz_divexact_ui(rden, q, g);
        }
    }
}

void fmpq_mul_ui(fmpq_t res, const fmpq_t op1, ulong c)
{
    _fmpq_mul_ui(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1), c);
}
