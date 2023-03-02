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

void check(ulong n, mp_limb_t ans)
{
    int ok, reasonable;
    mp_limb_t lo, hi;
    n_nth_prime_bounds(&lo, &hi, n);

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
    
    flint_printf("nth_prime_bounds....");
    fflush(stdout);

    for (n=6; n<7500 * FLINT_MIN(10, flint_test_multiplier()); n++)
    {
        check(n, n_nth_prime(n));
    }

    /* Some known large primes */
    check(UWORD(10), UWORD(29));
    check(UWORD(100), UWORD(541));
    check(UWORD(1000), UWORD(7919));
    check(UWORD(10000), UWORD(104729));
    check(UWORD(100000), UWORD(1299709));
    check(UWORD(1000000), UWORD(15485863));
    check(UWORD(10000000), UWORD(179424673));
    check(UWORD(100000000), UWORD(2038074743));
#if FLINT64
    check(UWORD(1000000000), UWORD(22801763489));
    check(UWORD(10000000000), UWORD(252097800623));
    check(UWORD(100000000000), UWORD(2760727302517));
    check(UWORD(1000000000000), UWORD(29996224275833));
    check(UWORD(10000000000000), UWORD(323780508946331));
    check(UWORD(100000000000000), UWORD(3475385758524527));
    check(UWORD(1000000000000000), UWORD(37124508045065437));
    check(UWORD(10000000000000000), UWORD(394906913903735329));
    check(UWORD(100000000000000000), UWORD(4185296581467695669));
#endif

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
