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
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"
#include "profiler.h"


int main(void)
{
    long i;
    ulong n;
    fmpz_t x, y, z;
    flint_rand_t state;

    printf("euler_phi....");
    fflush(stdout);

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(z);
    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        fmpz_set_ui(x, i);
        arith_euler_phi(y, x);
        arith_euler_phi(x, x);
        fmpz_set_ui(z, n_euler_phi(i));
        if (!fmpz_equal(x, y) || !fmpz_equal(x, z))
        {
            printf("FAIL: %ld\n", i);
            abort();
        }
    }

    /* Aliasing test */
    for (i = 0; i < 1000; i++)
    {
        fmpz_randtest(x, state, FLINT_BITS);
        fmpz_randtest(y, state, 5);
        fmpz_pow_ui(y, y, n_randtest(state) % 100);
        fmpz_mul(x, x, y);
        fmpz_set(z, x);
        arith_euler_phi(y, x);
        arith_euler_phi(x, x);
        if (!fmpz_equal(x, y))
        {
            printf("FAIL: ");
            fmpz_print(z);
            printf("\n");
            abort();
        }
    }

    /* Power of a single prime, phi(p^n) = (p-1) * p^(n-1) */
    n_compute_primes(100);
    for (i = 0; i < 100; i++)
    {
        n = (n_randtest(state) % 100) + 1;
        fmpz_set_ui(x, flint_primes[i]);
        fmpz_pow_ui(x, x, n);
        arith_euler_phi(x, x);
        fmpz_set_ui(y, flint_primes[i]);
        fmpz_pow_ui(y, y, n-1);
        fmpz_mul_ui(y, y, flint_primes[i]-1);
        if (!fmpz_equal(x, y))
        {
            printf("FAIL: %lu ^ %lu\n", flint_primes[i], n);
        }
    }

    /* Something nontrivial */
    fmpz_set_str(x, "10426024348053113487152988625265848110501553295256578345594388516660144", 10);
    fmpz_set_str(y, "2265085829098571747262267425315881590169106756213617459200000000000000", 10);
    arith_euler_phi(x, x);
    if (!fmpz_equal(x, y))
    {
        printf("FAIL: special test value\n");
        abort();
    }

    flint_randclear(state);

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(z);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
