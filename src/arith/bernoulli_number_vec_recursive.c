/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "arith.h"

static void
__ramanujan_even_common_denom(fmpz * num, fmpz * den, slong start, slong n)
{
    fmpz_t t, c, d, cden;
    slong j, k, m, mcase;
    int prodsize = 0;

    if (start >= n)
        return;

    fmpz_init(t);
    fmpz_init(c);
    fmpz_init(d);
    fmpz_init(cden);

    /* Common denominator */
    arith_primorial(cden, n + 1);

    start += start % 2;

    /* Convert initial values to common denominator */
    for (k = 0; k < start; k += 2)
    {
        fmpz_divexact(t, cden, den + k);
        fmpz_mul(num + k, num + k, t);
    }

    /* Ramanujan's recursive formula */
    for (m = start; m < n; m += 2)
    {
        mcase = m % 6;

        fmpz_mul_ui(num + m, cden, m + UWORD(3));
        fmpz_divexact_ui(num + m, num + m, UWORD(3));

        if (mcase == 4)
        {
            fmpz_neg(num + m, num + m);
            fmpz_divexact_ui(num + m, num + m, UWORD(2));
        }

        /* All factors are strictly smaller than m + 4; choose prodsize such
           that (m + 4)^prodsize fits in an slong. */
        {
#if FLINT64
            if      (m < WORD(1444))       prodsize = 6;
            else if (m < WORD(2097148))    prodsize = 3;
            else if (m < WORD(3037000495)) prodsize = 2;  /* not very likely... */
            else flint_throw(FLINT_ERROR, "(%s)", __func__);
#else
            if      (m < WORD(32))    prodsize = 6;
            else if (m < WORD(1286))  prodsize = 3;
            else if (m < WORD(46336)) prodsize = 2;
            else flint_throw(FLINT_ERROR, "(%s)", __func__);
#endif
        }

        /* c = t = binomial(m+3, m) */
        fmpz_set_ui(t, m + UWORD(1));
        fmpz_mul_ui(t, t, m + UWORD(2));
        fmpz_mul_ui(t, t, m + UWORD(3));
        fmpz_divexact_ui(t, t, UWORD(6));
        fmpz_set(c, t);

        for (j = 6; j <= m; j += 6)
        {
            slong r = m - j;

            /* c = binomial(m+3, m-j); */
            switch (prodsize)
            {
                case 2:
                fmpz_mul_ui(c, c, (r+6)*(r+5));
                fmpz_mul_ui(c, c, (r+4)*(r+3));
                fmpz_mul_ui(c, c, (r+2)*(r+1));
                fmpz_set_ui(d,    (j+0)*(j+3));
                fmpz_mul_ui(d, d, (j-2)*(j+2));
                fmpz_mul_ui(d, d, (j-1)*(j+1));
                fmpz_divexact(c, c, d);
                break;

                case 3:
                fmpz_mul_ui(c, c, (r+6)*(r+5)*(r+4));
                fmpz_mul_ui(c, c, (r+3)*(r+2)*(r+1));
                fmpz_set_ui(d,    (j+0)*(j+3)*(j-2));
                fmpz_mul_ui(d, d, (j+2)*(j-1)*(j+1));
                fmpz_divexact(c, c, d);
                break;

                case 6:
                fmpz_mul_ui(c, c,      (r+6)*(r+5)*(r+4)*(r+3)*(r+2)*(r+1));
                fmpz_divexact_ui(c, c, (j+0)*(j+3)*(j-2)*(j+2)*(j-1)*(j+1));
                break;
            }

            fmpz_submul(num + m, c, num + (m - j));
        }
        fmpz_divexact(num + m, num + m, t);
    }

    /* Convert to separate denominators */
    for (k = 0; k < n; k += 2)
    {
        arith_bernoulli_number_denom(den + k, k);
        fmpz_divexact(t, cden, den + k);
        fmpz_divexact(num + k, num + k, t);
    }

    fmpz_clear(t);
    fmpz_clear(c);
    fmpz_clear(d);
    fmpz_clear(cden);
}

void _arith_bernoulli_number_vec_recursive(fmpz * num, fmpz * den, slong n)
{
    slong i, start;
    fmpz_t t;
    fmpz_t d;

    fmpz_init(t);
    fmpz_init(d);

    start = FLINT_MIN(BERNOULLI_SMALL_NUMER_LIMIT, n);

    /* Initial values */
    for (i = 0; i < start; i += 2)
        _arith_bernoulli_number(num + i, den + i, i);

    __ramanujan_even_common_denom(num, den, start, n);

    /* Odd values */
    for (i = 1; i < n; i += 2)
        _arith_bernoulli_number(num + i, den + i, i);

    fmpz_clear(d);
    fmpz_clear(t);
}
