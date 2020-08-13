/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

#define BERNOULLI_DENOM_MAX_SMALL 178

#if FLINT64
#define __u32 unsigned int
#else
#define __u32 mp_limb_t
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

void arith_bernoulli_number_denom(fmpz_t den, ulong n)
{
    slong i;
    mp_limb_t p;
    const mp_limb_t * primes;

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
        n_prime_pi_bounds(&p, &p, n);
        primes = n_primes_arr_readonly(p + 2);

        fmpz_set_ui(den, UWORD(6));
        for (i = 2; i < n; i++)
        {
            p = primes[i];
            if (p - 1 > n)
                break;
            if (n % (p - 1) == 0)
                fmpz_mul_ui(den, den, p);
        }
    }
}
