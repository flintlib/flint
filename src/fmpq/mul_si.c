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
_fmpq_mul_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            slong r)
{
    if (r == 0 || fmpz_is_zero(p))
    {
        fmpz_zero(rnum);
        fmpz_one(rden);
    }
    else if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && r >= COEFF_MIN && r <= COEFF_MAX)
    {
        _fmpq_mul_small(rnum, rden, *p, *q, r, 1);
    }
    else if (r == 1)
    {
        fmpz_set(rnum, p);
        fmpz_set(rden, q);
    }
    else if (r == -WORD(1))
    {
        fmpz_neg(rnum, p);
        fmpz_set(rden, q);
    }
    else
    {
        ulong a, g;

        a = FLINT_ABS(r);
        g = _fmpz_gcd_ui(q, a);

        if (g == 1)
        {
            fmpz_set(rden, q);
            fmpz_mul_si(rnum, p, r);
        }
        else
        {
            /* not using fmpz_mul_si(...)  because of the special case g = -WORD_MIN */
            fmpz_mul_ui(rnum, p, a / g);
            if (r < 0)
                fmpz_neg(rnum, rnum);
            fmpz_divexact_ui(rden, q, g);
        }
    }
}

void fmpq_mul_si(fmpq_t res, const fmpq_t op1, slong c)
{
    _fmpq_mul_si(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1), c);
}
