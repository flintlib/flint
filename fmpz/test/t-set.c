/*
    Copyright (C) 2009 William Hart

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
#include "fmpz.h"

int
main(void)
{
    int i, result;

    FLINT_TEST_INIT(state);
    

    flint_printf("set....");
    fflush(stdout);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t c, d;
        flint_bitcnt_t bits;

        mpz_init(c);
        mpz_init(d);

        bits = n_randint(state, 200) + 1;

        _flint_rand_init_gmp(state);
        mpz_rrandomb(c, state->gmp_state, bits);

        if (n_randint(state, 2))
            mpz_neg(c, c);

        fmpz_init(a);
        fmpz_init(b);

        fmpz_set_mpz(a, c);
        fmpz_set(b, a);
        fmpz_get_mpz(d, b);

        result = (mpz_cmp(c, d) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("c = %Zd, d = %Zd\n", c, d);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(c);
        mpz_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
