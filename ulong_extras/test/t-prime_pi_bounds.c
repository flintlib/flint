/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

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

    printf("FAIL:\n");
    printf("n = %lu: %lu < %lu < %lu\n", n, lo, ans, hi);
    abort();
}

int main(void)
{
    int n;

    printf("prime_pi_bounds....");
    fflush(stdout);

    for (n=17; n<10000 * FLINT_MIN(10, flint_test_multiplier()); n++)
    {
        check(n, n_prime_pi(n));
    }

    check(10UL, 4UL);
    check(100UL, 25UL);
    check(1000UL, 168UL);
    check(10000UL, 1229UL);
    check(100000UL, 9592UL);
    check(1000000UL, 78498UL);
    check(10000000UL, 664579UL);
    check(100000000UL, 5761455UL);
    check(1000000000UL, 50847534UL);
#if FLINT64
    check(10000000000UL, 455052511UL);
    check(100000000000UL, 4118054813UL);
    check(1000000000000UL, 37607912018UL);
    check(10000000000000UL, 346065536839UL);
    check(100000000000000UL, 3204941750802UL);
    check(1000000000000000UL, 29844570422669UL);
    check(10000000000000000UL, 279238341033925UL);
    check(100000000000000000UL, 2623557157654233UL);
#endif

    printf("PASS\n");
    return 0;
}
