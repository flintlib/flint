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
    

    flint_printf("get/set_mpz....");
    fflush(stdout);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t b, c;
        flint_bitcnt_t bits;

        mpz_init(b);
        mpz_init(c);

        bits = n_randint(state, 200) + 1;

        _flint_rand_init_gmp(state);
        mpz_rrandomb(b, state->gmp_state, bits);

        if (n_randint(state, 2))
            mpz_neg(b, b);

        fmpz_init(a);

        fmpz_set_mpz(a, b);
        fmpz_get_mpz(c, a);

        result = (mpz_cmp(b, c) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("b = %Zd, c = %Zd\n", b, c);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);

        mpz_clear(b);
        mpz_clear(c);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
