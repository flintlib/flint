/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2016 William Hart

    This file is part of FLINT.
 
   FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

int
_fmpq_cmp_ui(const fmpz_t p, const fmpz_t q, ulong c)
{
    int s1, s2, res;
    flint_bitcnt_t bp, bq, bc;
    fmpz_t u;

    if (fmpz_is_one(q))
        return fmpz_cmp_ui(p, c);

    if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && c <= WORD_MAX)
    {
        ulong a1, a0, b1, b0;

        a0 = *p;
        a1 = FLINT_SIGN_EXT(a0);
        smul_ppmm(b1, b0, *q, c);
        sub_ddmmss(a1, a0, a1, a0, b1, b0);

        if ((slong) a1 < 0)
            return -1;
        if ((slong) a1 > 0)
            return 1;
        return a0 != 0;
    }

    s1 = fmpz_sgn(p);
    s2 = (c > 0);

    if (s1 != s2)
        return s1 < s2 ? -1 : 1;

    if (s1 == 0)
        return -s2;

    if (s2 == 0)
        return -s1;

    bp = fmpz_bits(p);
    bq = fmpz_bits(q);
    bc = FLINT_BIT_COUNT(c);

    if (bp + 2 < bc + bq)
        return -s1;

    if (bp > bc + bq)
        return s1;

    fmpz_init(u);

    fmpz_mul_ui(u, q, c);

    res = fmpz_cmp(p, u);

    fmpz_clear(u);

    return res;
}

int
fmpq_cmp_ui(const fmpq_t x, ulong c)
{
    return _fmpq_cmp_ui(fmpq_numref(x), fmpq_denref(x), c);
}

