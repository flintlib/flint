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
_fmpq_next_calkin_wilf(fmpz_t rnum, fmpz_t rden,
    const fmpz_t num, const fmpz_t den)
{
    fmpz n, d;

    n = *num;
    d = *den;

    if (!COEFF_IS_MPZ(n) && !COEFF_IS_MPZ(d))
    {
        /* This does not overflow, as the larger part at most doubles */
        fmpz_set_ui(rnum, d);
        fmpz_set_ui(rden, d*(n / d) + d - (n % d));
    }
    else
    {
        fmpz_t q, r, t;

        fmpz_init(q);
        fmpz_init(r);
        fmpz_init(t);

        fmpz_fdiv_qr(q, r, num, den);
        fmpz_set(rnum, den);
        fmpz_mul(t, q, den);
        fmpz_add(rden, t, den);
        fmpz_sub(rden, rden, r);

        fmpz_clear(q);
        fmpz_clear(r);
        fmpz_clear(t);
    }
}

void
fmpq_next_calkin_wilf(fmpq_t res, const fmpq_t x)
{
    _fmpq_next_calkin_wilf(fmpq_numref(res), fmpq_denref(res),
        fmpq_numref(x), fmpq_denref(x));
}
