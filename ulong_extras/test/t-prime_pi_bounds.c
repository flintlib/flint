/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"

void check(mp_limb_t n, ulong ans)
{
    int ok, reasonable;
    ulong lo, hi;
    n_prime_pi_bounds(&lo, &hi, n);

    ok = lo <= ans && ans <= hi;
    reasonable = (n < 1000) || (ans/2 < lo && hi < ans*2);

    if (ok && reasonable)
        return;

    flint_printf("FAIL:\n");
    flint_printf("n = %wu: %wu < %wu < %wu\n", n, lo, ans, hi);
    fflush(stdout);
    flint_abort();
}

int main(void)
{
    int n;

    FLINT_TEST_INIT(state);
    flint_printf("prime_pi_bounds....");
    fflush(stdout);

    for (n=17; n<10000 * FLINT_MIN(10, flint_test_multiplier()); n++)
    {
        check(n, n_prime_pi(n));
    }

    check(UWORD(10), UWORD(4));
    check(UWORD(100), UWORD(25));
    check(UWORD(1000), UWORD(168));
    check(UWORD(10000), UWORD(1229));
    check(UWORD(100000), UWORD(9592));
    check(UWORD(1000000), UWORD(78498));
    check(UWORD(10000000), UWORD(664579));
    check(UWORD(100000000), UWORD(5761455));
    check(UWORD(1000000000), UWORD(50847534));
#if FLINT64
    check(UWORD(10000000000), UWORD(455052511));
    check(UWORD(100000000000), UWORD(4118054813));
    check(UWORD(1000000000000), UWORD(37607912018));
    check(UWORD(10000000000000), UWORD(346065536839));
    check(UWORD(100000000000000), UWORD(3204941750802));
    check(UWORD(1000000000000000), UWORD(29844570422669));
    check(UWORD(10000000000000000), UWORD(279238341033925));
    check(UWORD(100000000000000000), UWORD(2623557157654233));
#endif

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
