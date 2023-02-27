/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q)
{
    if (fmpz_is_zero(p))
    {
        fmpq_zero(res);
    }
    else if (fmpz_is_pm1(q) || fmpz_is_pm1(p))
    {
        if (fmpz_sgn(q) < 0)
        {
            fmpz_neg(fmpq_numref(res), p);
            fmpz_neg(fmpq_denref(res), q);
        }
        else
        {
            fmpz_set(fmpq_numref(res), p);
            fmpz_set(fmpq_denref(res), q);
        }
    }
    else
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_gcd(t, p, q);

        if (fmpz_is_one(t))
        {
            fmpz_set(fmpq_numref(res), p);
            fmpz_set(fmpq_denref(res), q);
        }
        else
        {
            fmpz_divexact(fmpq_numref(res), p, t);
            fmpz_divexact(fmpq_denref(res), q, t);
        }

        if (fmpz_sgn(fmpq_denref(res)) < 0)
        {
            fmpz_neg(fmpq_numref(res), fmpq_numref(res));
            fmpz_neg(fmpq_denref(res), fmpq_denref(res));
        }

        fmpz_clear(t);
    }
}
