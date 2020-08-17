/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

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

    flint_printf("powm....");
    fflush(stdout);

    

    /* Compare with MPIR */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, m;
        fmpz_t x;
        mpz_t y;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(x);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(m);
        mpz_init(y);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(c, state, 200);
        fmpz_abs(c, c);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(m, c);
        fmpz_randtest_unsigned(x, state, 20);
        fmpz_get_mpz(y, x);

        fmpz_powm(b, a, x, c);
        mpz_powm(e, d, y, m);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(e, f) == 0);
        if (!result)
        {
            flint_printf("FAIL (cmp f with MPIR e := d^y mod m):\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, y = %Zd, m = %Zd\n", d, e, f, y, m);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(x);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(m);
        mpz_clear(y);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        fmpz_t n;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(n);

        fmpz_randtest(b, state, 200);
        fmpz_randtest_not_zero(c, state, 200);
        fmpz_abs(c, c);
        fmpz_randtest_unsigned(n, state, 20);

        fmpz_powm(a, b, n, c);
        fmpz_powm(b, b, n, c);

        result = (fmpz_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL (alias a and b):\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("n = "), fmpz_print(n), flint_printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(n);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        fmpz_t n;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(n);

        fmpz_randtest(b, state, 200);
        fmpz_randtest_not_zero(c, state, 200);
        fmpz_abs(c, c);
        fmpz_randtest_unsigned(n, state, 20);

        fmpz_powm(a, b, n, c);
        fmpz_powm(c, b, n, c);

        result = (fmpz_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (alias a and c):\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("n = "), fmpz_print(n), flint_printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(n);
    }

    /* Check aliasing of a and {b, c} */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        fmpz_t n;

        fmpz_init(a);
        fmpz_init(c);
        fmpz_init(n);

        fmpz_randtest_not_zero(c, state, 200);
        fmpz_abs(c, c);
        fmpz_randtest_unsigned(n, state, 20);

        fmpz_powm(a, c, n, c);
        fmpz_powm(c, c, n, c);

        result = (fmpz_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (alias a and b,c):\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("n = "), fmpz_print(n), flint_printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        fmpz_clear(n);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

