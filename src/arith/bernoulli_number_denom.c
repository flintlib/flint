/*
    Copyright (C) 2011, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "arith.h"

#define BERNOULLI_DENOM_MAX_SMALL 178

#if FLINT64
#define __u32 unsigned int
#else
#define __u32 ulong
#endif

static const __u32 __bernoulli_denom_small[] =
{
    UWORD(1), UWORD(6), UWORD(30), UWORD(42), UWORD(30), UWORD(66), UWORD(2730), UWORD(6), UWORD(510), UWORD(798), UWORD(330),
    UWORD(138), UWORD(2730), UWORD(6), UWORD(870), UWORD(14322), UWORD(510), UWORD(6), UWORD(1919190), UWORD(6), UWORD(13530),
    UWORD(1806), UWORD(690), UWORD(282), UWORD(46410), UWORD(66), UWORD(1590), UWORD(798), UWORD(870), UWORD(354),
    UWORD(56786730), UWORD(6), UWORD(510), UWORD(64722), UWORD(30), UWORD(4686), UWORD(140100870), UWORD(6), UWORD(30),
    UWORD(3318), UWORD(230010), UWORD(498), UWORD(3404310), UWORD(6), UWORD(61410), UWORD(272118), UWORD(1410), UWORD(6),
    UWORD(4501770), UWORD(6), UWORD(33330), UWORD(4326), UWORD(1590), UWORD(642), UWORD(209191710), UWORD(1518),
    UWORD(1671270), UWORD(42), UWORD(1770), UWORD(6), UWORD(2328255930), UWORD(6), UWORD(30), UWORD(4357878), UWORD(510),
    UWORD(8646), UWORD(4206930), UWORD(6), UWORD(4110), UWORD(274386), UWORD(679470), UWORD(6), UWORD(2381714790),
    UWORD(6), UWORD(4470), UWORD(2162622), UWORD(30), UWORD(138), UWORD(1794590070), UWORD(6), UWORD(230010),
    UWORD(130074), UWORD(2490), UWORD(1002), UWORD(3404310), UWORD(66), UWORD(5190), UWORD(2478), UWORD(1043970),
    UWORD(1074),
};

/* For even n > 178: the product of the primes d + 1 over the divisors d
   of n (von Staudt-Clausen), the divisors enumerated from the
   factorization of n -- a few hundred primality tests at most for
   n < 2^32, instead of a scan of the primes up to n + 1. */
void arith_bernoulli_number_denom(fmpz_t den, ulong n)
{
    if (n % 2 == 1)
    {
        fmpz_set_ui(den, 1 + (n == 1));
    }
    else if (n <= BERNOULLI_DENOM_MAX_SMALL)
    {
        fmpz_set_ui(den, __bernoulli_denom_small[n / 2]);
    }
    else
    {
        n_factor_t fac;
        ulong * divs;
        slong i, j, e, nd, m, ndiv;

        n_factor_init(&fac);
        n_factor(&fac, n, 1);

        ndiv = 1;
        for (i = 0; i < fac.num; i++)
            ndiv *= fac.exp[i] + 1;
        divs = flint_malloc(ndiv * sizeof(ulong));

        divs[0] = 1;
        nd = 1;
        for (i = 0; i < fac.num; i++)
        {
            ulong p = fac.p[i], pe = 1;

            m = nd;
            for (e = 1; e <= fac.exp[i]; e++)
            {
                pe *= p;
                for (j = 0; j < m; j++)
                    divs[nd++] = divs[j] * pe;
            }
        }

        fmpz_one(den);
        for (i = 0; i < nd; i++)
            if (n_is_prime(divs[i] + 1))
                fmpz_mul_ui(den, den, divs[i] + 1);

        flint_free(divs);
    }
}
