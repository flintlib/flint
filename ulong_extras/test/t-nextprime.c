/*
    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    mp_limb_t n;
    mp_limb_t res1, res2;
    slong rep;
    mpz_t mpz_n;
    FLINT_TEST_INIT(state);
    
    flint_printf("nextprime....");
    fflush(stdout);

    if (n_nextprime(0, 0) != 2)
    {
        flint_printf("FAIL: expected n_nextprime(0) = 2");
        abort();
    }

    if (n_nextprime(UWORD_MAX_PRIME - 1, 0) != UWORD_MAX_PRIME)
    {
        flint_printf("FAIL: expected n_nextprime(UWORD_MAX_PRIME-1) = UWORD_MAX_PRIME");
        abort();
    }

    mpz_init(mpz_n);

    for (rep = 0; rep < 10000 * flint_test_multiplier(); rep++)
    {
        ulong bits = n_randint(state, FLINT_D_BITS-1)+1;
        n = n_randtest(state) % ((UWORD(1)<<bits) - UWORD(1)) + 1; 
        
        flint_mpz_set_ui(mpz_n, n);

        mpz_nextprime(mpz_n, mpz_n);
        n = n_nextprime(n, 0);

        res1 = n;
        res2 = flint_mpz_get_ui(mpz_n);

        if (res1 != res2)
        {
            flint_printf("FAIL:\n");
            flint_printf("%wu, %wu\n", res1, res2); 
            abort();
        }
    }

    for (rep = 0; rep < 10000; rep++)
    {
        n = (UWORD(1) << (FLINT_BITS-1)) + rep;
        flint_mpz_set_ui(mpz_n, n);

        mpz_nextprime(mpz_n, mpz_n);
        n = n_nextprime(n, 0);

        res1 = n;
        res2 = flint_mpz_get_ui(mpz_n);

        if (res1 != res2)
        {
            flint_printf("FAIL:\n");
            flint_printf("%wu, %wu\n", res1, res2); 
            abort();
        }
    }

    mpz_clear(mpz_n); 

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
