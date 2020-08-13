/*
    Copyright (C) 2010 Sebastian Pancratz

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
#include "long_extras.h"
#include "fmpz.h"

int main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    
    flint_printf("mpz_init_set_readonly....");
    fflush(stdout);
    
    

    /* Create some small fmpz integers, clear the mpz_t */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t z;

        *f = z_randint(state, COEFF_MAX + 1);

        flint_mpz_init_set_readonly(z, f);
        flint_mpz_clear_readonly(z);
    }

    /* Create some large fmpz integers, do not clear the mpz_t */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t z;

        fmpz_init(f);
        fmpz_randtest(f, state, 2 * FLINT_BITS);

        if (COEFF_IS_MPZ(*f))
        {
            flint_mpz_init_set_readonly(z, f);
        }

        fmpz_clear(f);
    }

    /* Create some more fmpz integers */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, g;
        mpz_t z;

        fmpz_init(f);
        fmpz_init(g);
        fmpz_randtest(f, state, 2 * FLINT_BITS);

        flint_mpz_init_set_readonly(z, f);
        fmpz_set_mpz(g, z);

        if (!fmpz_equal(f, g))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            gmp_printf("z = %Zd\n", z);
        }

        flint_mpz_clear_readonly(z);
	fmpz_clear(f);
        fmpz_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
