/*
    Copyright (C) 2012, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

int
_fmpq_cmp(const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s)
{
    int s1, s2, res;
    flint_bitcnt_t bp, bq, br, bs;
    fmpz_t t, u;

    if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && !COEFF_IS_MPZ(*r) && !COEFF_IS_MPZ(*s))
    {
        slong a1, a0, b1, b0;

        smul_ppmm(a1, a0, *p, *s);
        smul_ppmm(b1, b0, *q, *r);
        sub_ddmmss(a1, a0, a1, a0, b1, b0);

        if (a1 < 0)
            return -1;
        if (a1 > 0)
            return 1;
        return a0 != 0;
    }

    if (fmpz_equal(q, s))
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
    bs = fmpz_bits(s);

    if (bp + bs + 1 < br + bq)
        return -s1;

    if (bp + bs > br + bq + 1)
        return s1;

    if (fmpz_is_one(q))
    {
        fmpz_init(t);
        fmpz_mul(t, p, s);
        res = fmpz_cmp(t, r);
        fmpz_clear(t);
    }
    else if (fmpz_is_one(s))
    {
        fmpz_init(u);
        fmpz_mul(u, q, r);

        res = fmpz_cmp(p, u);

        fmpz_clear(u);
    }
    else
    {
        fmpz_init(t);
        fmpz_init(u);

        fmpz_mul(t, p, s);
        fmpz_mul(u, q, r);

        res = fmpz_cmp(t, u);

        fmpz_clear(t);
        fmpz_clear(u);
    }

    return res;
}

int
fmpq_cmp(const fmpq_t x, const fmpq_t y)
{
    return _fmpq_cmp(fmpq_numref(x), fmpq_denref(x),
                     fmpq_numref(y), fmpq_denref(y));
}

