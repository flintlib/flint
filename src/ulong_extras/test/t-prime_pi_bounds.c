/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

void check_prime_pi_bound(mp_limb_t n, ulong ans)
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

TEST_FUNCTION_START(n_prime_pi_bounds, state)
{
    int n;

    for (n=17; n<10000 * FLINT_MIN(10, flint_test_multiplier()); n++)
    {
        check_prime_pi_bound(n, n_prime_pi(n));
    }

    check_prime_pi_bound(UWORD(10), UWORD(4));
    check_prime_pi_bound(UWORD(100), UWORD(25));
    check_prime_pi_bound(UWORD(1000), UWORD(168));
    check_prime_pi_bound(UWORD(10000), UWORD(1229));
    check_prime_pi_bound(UWORD(100000), UWORD(9592));
    check_prime_pi_bound(UWORD(1000000), UWORD(78498));
    check_prime_pi_bound(UWORD(10000000), UWORD(664579));
    check_prime_pi_bound(UWORD(100000000), UWORD(5761455));
    check_prime_pi_bound(UWORD(1000000000), UWORD(50847534));
#if FLINT64
    check_prime_pi_bound(UWORD(10000000000), UWORD(455052511));
    check_prime_pi_bound(UWORD(100000000000), UWORD(4118054813));
    check_prime_pi_bound(UWORD(1000000000000), UWORD(37607912018));
    check_prime_pi_bound(UWORD(10000000000000), UWORD(346065536839));
    check_prime_pi_bound(UWORD(100000000000000), UWORD(3204941750802));
    check_prime_pi_bound(UWORD(1000000000000000), UWORD(29844570422669));
    check_prime_pi_bound(UWORD(10000000000000000), UWORD(279238341033925));
    check_prime_pi_bound(UWORD(100000000000000000), UWORD(2623557157654233));
#endif

    TEST_FUNCTION_END(state);
}
