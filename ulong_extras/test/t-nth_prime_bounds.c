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

void check(ulong n, mp_limb_t ans)
{
    int ok, reasonable;
    mp_limb_t lo, hi;
    n_nth_prime_bounds(&lo, &hi, n);

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

    printf("nth_prime_bounds....");
    fflush(stdout);

    for (n=6; n<7500 * FLINT_MIN(10, flint_test_multiplier()); n++)
    {
        check(n, n_nth_prime(n));
    }

    /* Some known large primes */
    check(10UL, 29UL);
    check(100UL, 541UL);
    check(1000UL, 7919UL);
    check(10000UL, 104729UL);
    check(100000UL, 1299709UL);
    check(1000000UL, 15485863UL);
    check(10000000UL, 179424673UL);
    check(100000000UL, 2038074743UL);
#if FLINT64
    check(1000000000UL, 22801763489UL);
    check(10000000000UL, 252097800623UL);
    check(100000000000UL, 2760727302517UL);
    check(1000000000000UL, 29996224275833UL);
    check(10000000000000UL, 323780508946331UL);
    check(100000000000000UL, 3475385758524527UL);
    check(1000000000000000UL, 37124508045065437UL);
    check(10000000000000000UL, 394906913903735329UL);
    check(100000000000000000UL, 4185296581467695669UL);
#endif

    printf("PASS\n");
    return 0;
}
