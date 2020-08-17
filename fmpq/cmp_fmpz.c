/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"


int _fmpq_cmp_fmpz(const fmpz_t p, const fmpz_t q, const fmpz_t r)
{
    int s1, s2, res;
    flint_bitcnt_t bp, bq, br, bs;
    fmpz_t u;

    if (fmpz_is_one(q))
        return fmpz_cmp(p, r);

    s1 = fmpz_sgn(p);
    s2 = fmpz_sgn(r);

    if (s1 != s2)
        return s1 < s2 ? -1 : 1;

    if (s1 == 0)
        return -s2;

    if (s2 == 0)
        return -s1;

    bp = fmpz_bits(p);
    bq = fmpz_bits(q);
    br = fmpz_bits(r);
    bs = 1;

    if (bp + bs + 1 < br + bq)
        return -s1;

    if (bp + bs > br + bq + 1)
        return s1;

    fmpz_init(u);

    fmpz_mul(u, q, r);

    res = fmpz_cmp(p, u);

    fmpz_clear(u);

    return res;
}


int fmpq_cmp_fmpz(const fmpq_t x, const fmpz_t y)
{
    return _fmpq_cmp_fmpz(fmpq_numref(x), fmpq_denref(x), y);
}

