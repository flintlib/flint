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

/* Use the definiton of GMP versions >= 6.0 */
int
mpz_invert2(mpz_t a, const mpz_t b, const mpz_t c)
{
    if (mpz_cmpabs_ui(c, 1) == 0)
    {
        mpz_set_ui(a, 0);
        return 1;
    }
    else
        return mpz_invert(a, b, c);
}

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("invmod....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, g;
        int r1, r2;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        r1 = fmpz_invmod(c, a, b);
        r2 = mpz_invert2(f, d, e);

        fmpz_get_mpz(g, c);

        result = (r1 != 0 && r2 != 0 && (mpz_cmp(f, g) == 0)) || (r1 == 0 && r2 == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, r1 = %d, r2 = %d\n", d,
                 e, f, g, r1, r2);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        mpz_t d, f, g;
        int r1, r2;

        fmpz_init(a);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest_not_zero(a, state, 200);

        fmpz_get_mpz(d, a);

        r1 = fmpz_invmod(c, a, a);
        r2 = mpz_invert2(f, d, d);

        fmpz_get_mpz(g, c);

        result = (r1 != 0 && r2 != 0 && (mpz_cmp(f, g) == 0)) || (r1 == 0 && r2 == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, f = %Zd, g = %Zd, r1 = %d, r2 = %d\n", d, f,
                       g, r1, r2);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Test aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f, g;
        int r1, r2;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        r1 = fmpz_invmod(a, a, b);
        r2 = mpz_invert2(f, d, e);

        fmpz_get_mpz(g, a);

        result = (r1 != 0 && r2 != 0 && (mpz_cmp(f, g) == 0)) || (r1 == 0 && r2 == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd\n", d, e, f, g);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Test aliasing of b and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f, g;
        int r1, r2;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        r1 = fmpz_invmod(b, a, b);
        r2 = mpz_invert2(f, d, e);

        fmpz_get_mpz(g, b);

        result = (r1 != 0 && r2 != 0 && (mpz_cmp(f, g) == 0)) || (r1 == 0 && r2 == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd\n", d, e, f, g);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
