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

    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    mp_limb_t n;
    mp_limb_t res1, res2;
    long i, rep;
    mpz_t mpz_n;
    flint_rand_t state;
    
    printf("nextprime....");
    fflush(stdout);

    flint_randinit(state);

    if (n_nextprime(0, 0) != 2)
    {
        printf("FAIL: expected n_nextprime(0) = 2");
        abort();
    }

    if (n_nextprime(ULONG_MAX_PRIME - 1, 0) != ULONG_MAX_PRIME)
    {
        printf("FAIL: expected n_nextprime(ULONG_MAX_PRIME-1) = ULONG_MAX_PRIME");
        abort();
    }

    mpz_init(mpz_n);

    for (rep = 0; rep < 10000 * flint_test_multiplier(); rep++)
    {
        unsigned long bits = n_randint(state, FLINT_D_BITS-1)+1;
        n = n_randtest(state) % ((1UL<<bits) - 1UL) + 1; 
        mpz_set_ui(mpz_n, n);

        for (i = 0; i < 1; i++)
        {
            mpz_nextprime(mpz_n, mpz_n);
            n = n_nextprime(n, 0);
        }

        res1 = n;
        res2 = mpz_get_ui(mpz_n);

        if (res1 != res2)
        {
            printf("FAIL:\n");
            printf("%lu, %lu\n", res1, res2); 
            abort();
        }
    }

    mpz_clear(mpz_n); 

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
