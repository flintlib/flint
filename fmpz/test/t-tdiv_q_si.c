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
#include "fmpz.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("tdiv_q_si....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        slong b;
        fmpz_t a, c;
        mpz_t d, e, f, g;

        fmpz_init(a);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        b = z_randtest_not_zero(state);

        fmpz_get_mpz(d, a);
        flint_mpz_set_si(e, b);

        fmpz_tdiv_q_si(c, a, b);
        mpz_tdiv_q(f, d, e);

        fmpz_get_mpz(g, c);

        result = (mpz_cmp(f, g) == 0);
        if (!result)
        {
            flint_printf("FAIL (1):\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd\n", d, e, f, g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Test aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        slong b;
        fmpz_t a;
        mpz_t d, e, f, g;

        fmpz_init(a);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        b = z_randtest_not_zero(state);

        fmpz_get_mpz(d, a);
        flint_mpz_set_si(e, b);

        fmpz_tdiv_q_si(a, a, b);
        mpz_tdiv_q(f, d, e);

        fmpz_get_mpz(g, a);

        result = (mpz_cmp(f, g) == 0);
        if (!result)
        {
            flint_printf("FAIL (2):\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd\n", d, e, f, g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
