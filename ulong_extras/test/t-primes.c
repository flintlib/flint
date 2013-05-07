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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    flint_rand_t state;
    long n;

    printf("primes....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    /* compare with n_nextprime */
    {
        n_primes_t iter;
        long i;
        mp_limb_t p, q;

        n_primes_init(iter);
        q = 0;
        for (i = 0; i < 200000; i++)
        {
            p = n_primes_next(iter);
            q = n_nextprime(q, 0);

            if (p != q)
            {
                printf("FAIL\n");
                printf("i = %lu, p = %lu, q = %lu\n", i, p, q);
                abort();
            }
        }

        n_primes_clear(iter);
    }

    /* count primes */
    for (n = 0; n < 10; n++)
    {
        n_primes_t iter;
        mp_limb_t s, p, r;

        const unsigned int primepi[10] = {
            0, 4, 25, 168, 1229, 9592, 78498, 664579, 5761455, 50847534
        };

        r = n_pow(10, n);

        n_primes_init(iter);
        s = 0;
        while ((p = n_primes_next(iter)) <= r)
            s++;

        if (s != primepi[n])
        {
            printf("FAIL\n");
            printf("pi(10^%ld) = %u, computed = %lu\n", n, primepi[n], s);
            abort();
        }

        n_primes_clear(iter);
    }


    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
