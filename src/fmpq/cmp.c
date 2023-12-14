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
        ulong a1, a0, b1, b0;

        smul_ppmm(a1, a0, *p, *s);
        smul_ppmm(b1, b0, *q, *r);
        sub_ddmmss(a1, a0, a1, a0, b1, b0);

        if ((slong) a1 < 0)
            return -1;
        if ((slong) a1 > 0)
            return 1;
        return a0 != 0;
    }

    if (fmpz_equal(q, s))
        return fmpz_cmp(p, r);

    s1 = fmpz_sgn(p);
    s2 = fmpz_sgn(r);

    if (s1 != s2)
        return s1 < s2 ? -1 : 1;

    /* NOTE: If `p' and `r' was zero, then we stepped into the first
     * if-statement. Else, if `p' or `r' was zero, then s1 != s2. Hence, at this
     * stage `p' and `r' has to be non-zero. */

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

int _fmpq_cmp_fmpz(const fmpz_t p, const fmpz_t q, const fmpz_t r)
{
    fmpz one = 1;
    return _fmpq_cmp(p, q, r, &one);
}


int fmpq_cmp_fmpz(const fmpq_t x, const fmpz_t y)
{
    return _fmpq_cmp_fmpz(fmpq_numref(x), fmpq_denref(x), y);
}

int
_fmpq_cmp_si(const fmpz_t p, const fmpz_t q, slong c)
{
    int s1, s2, res;
    flint_bitcnt_t bp, bq, bc;
    slong d;
    fmpz_t u;

    if (fmpz_is_one(q))
        return fmpz_cmp_si(p, c);

    if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q))
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
    s2 = c > 0 ? 1 : (c < 0 ? -1 : 0);

    if (s1 != s2)
        return s1 < s2 ? -1 : 1;

    if (s1 == 0)
        return 0;

    bp = fmpz_bits(p);
    bq = fmpz_bits(q);

    d = -c;

    if (c != d) /* check for SLONG_MIN */
    {
        d = c < 0 ? -c : c;

        bc = FLINT_BIT_COUNT(d);

        if (bp + 2 < bc + bq)
            return -s1;

        if (bp > bc + bq)
            return s1;
    }

    fmpz_init(u);

    fmpz_mul_si(u, q, c);

    res = fmpz_cmp(p, u);

    fmpz_clear(u);

    return res;
}

int
fmpq_cmp_si(const fmpq_t x, slong c)
{
    return _fmpq_cmp_si(fmpq_numref(x), fmpq_denref(x), c);
}

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

    /* NOTE: If `p' and `r' was zero, then we stepped into the first
     * if-statement. Else, if `p' or `r' was zero, then s1 != s2. Hence, at this
     * stage `p' and `r' has to be non-zero. */

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
