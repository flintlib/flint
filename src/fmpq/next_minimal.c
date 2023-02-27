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
_fmpq_next_minimal(fmpz_t rnum, fmpz_t rden,
                        const fmpz_t num, const fmpz_t den)
{
    fmpz p, q;

    p = *num;
    q = *den;

    if (!COEFF_IS_MPZ(p) && !COEFF_IS_MPZ(q))
    {
        if (p < q && p)
        {
            fmpz_set_ui(rnum, q);
            fmpz_set_ui(rden, p);
            return;
        }

        while (q < p)
        {
            q++;
            if (n_gcd(p, q) == 1)
            {
                fmpz_set_ui(rnum, q);
                fmpz_set_ui(rden, p);
                return;
            }
        }

        fmpz_one(rnum);
        fmpz_set_ui(rden, p + 1);
    }
    else
    {
        fmpz_t t;

        if (fmpz_cmp(num, den) < 0)
        {
            fmpz_set(rnum, num);
            fmpz_set(rden, den);
            fmpz_swap(rnum, rden);
            return;
        }

        fmpz_init(t);
        fmpz_set(rnum, num);
        fmpz_set(rden, den);

        while (fmpz_cmp(rden, rnum) < 0)
        {
            fmpz_add_ui(rden, rden, 1);
            fmpz_gcd(t, rden, rnum);
            if (fmpz_is_one(t))
            {
                fmpz_swap(rnum, rden);
                fmpz_clear(t);
                return;
            }
        }

        fmpz_add_ui(rden, rden, 1);
        fmpz_one(rnum);
        fmpz_clear(t);
    }
}

void
fmpq_next_minimal(fmpq_t res, const fmpq_t x)
{
    _fmpq_next_minimal(fmpq_numref(res), fmpq_denref(res),
        fmpq_numref(x), fmpq_denref(x));
}
