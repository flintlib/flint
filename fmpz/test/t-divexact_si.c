/*
    Copyright (C) 2009 William Hart
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

    flint_printf("divexact_si....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        mpz_t e, f, g;
        slong n;

        fmpz_init(a);
        fmpz_init(c);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        n = z_randtest_not_zero(state);
        fmpz_mul_si(c, a, n);

        fmpz_get_mpz(e, c);

        fmpz_divexact_si(a, c, n);
        flint_mpz_divexact_ui(f, e, FLINT_ABS(n));
        if (n < 0)
            mpz_neg(f, f);

        fmpz_get_mpz(g, a);

        result = (mpz_cmp(f, g) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("n = %Md, e = %Zd, f = %Zd, g = %Zd\n", n, e, f, g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Test aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        mpz_t d, f, g;
        slong n;

        fmpz_init(a);
        fmpz_init(c);
        mpz_init(d);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        n = z_randtest_not_zero(state);
        fmpz_mul_si(c, a, n);

        fmpz_get_mpz(d, c);

        fmpz_divexact_si(c, c, n);
        flint_mpz_divexact_ui(f, d, FLINT_ABS(n));
        if (n < 0)
            mpz_neg(f, f);

        fmpz_get_mpz(g, c);

        result = (mpz_cmp(f, g) == 0);

        if (!result)
        {
            flint_printf("FAIL;\n");
            gmp_printf("d = %Zd, n = %Md, f = %Zd, g = %Zd\n", d, n, f, g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        mpz_clear(d);
        mpz_clear(f);
        mpz_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
