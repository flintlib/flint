/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_rootrem, state)
{
   int i, result;
   mp_limb_t upper_limit;

#if FLINT64
   upper_limit = 2642245;
#else
   upper_limit = 1625;
#endif

    /* random n and root */

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        mp_limb_t a, b, c, d, val, j;
        mpz_t e, f, g, h;

        mpz_init(e);
        mpz_init(f);
        mpz_init(g);
        mpz_init(h);

        c = n_randint(state, 0);    /*number */
        flint_mpz_set_ui(g, c);

        d = n_randint(state, 0);   /*root */
        flint_mpz_set_ui(h, d);

        a = n_rootrem(&b, c, d);

        mpz_rootrem(e, f, g, flint_mpz_get_ui(h));

        val = flint_mpz_get_ui(e);
        j = flint_mpz_get_ui(f);

        result = ((a == val) && (b == j));

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("Passed Parameters : n = %wu root = %wu", c, d);
            flint_printf("Answer generated : base = %wu remainder = %wu", a, b);
            flint_printf("Expected answer : base = %wu remainder = %wu", val, j);
            fflush(stdout);
            flint_abort();
        }
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
        mpz_clear(h);
    }

    /* n of type a^b */

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        mp_limb_t a, b, c, d, max_pow, base;

        base = n_randint(state, upper_limit - 2) + 2;     /* base form 2 to 2642245*/
        max_pow = n_flog(UWORD_MAX, base);
        d = n_randint(state, max_pow);       /* root */
        if (!d)
            d+=1;

        c = n_pow(base, d);                  /* number */
        a = n_rootrem(&b, c, d);
        result = ((a == base) && (b == 0));

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("Passed Parameters : n = %wu root = %wu", c, d);
            printf("\n");
            flint_printf("Answer generated : base = %wu remainder = %wu", a, b);
            printf("\n");
            flint_printf("Expected answer : base = %wu remainder = 0", base);
            fflush(stdout);
            flint_abort();
        }
    }

    /* n of type a^b + 1 */

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        mp_limb_t a, b, c, d, max_pow, base;

        base = n_randint(state, upper_limit - 2) + 2;     /* base between 2 to 2642245*/
        max_pow = n_flog(UWORD_MAX, base);
        d = n_randint(state, max_pow);
        if (d < 2)                                /* root between 2 to max_pow */
            d = 2;

        c = n_pow(base, d) + 1;                   /* number */
        a = n_rootrem(&b, c, d);
        result = ((a == base) && (b == 1));

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("Passed Parameters : n = %wu root = %wu", c, d);
            printf("\n");
            flint_printf("Answer generated : base = %wu remainder = %wu", a, b);
            printf("\n");
            flint_printf("Expected answer : base = %wu remainder = 1", base);
            fflush(stdout);
            flint_abort();
        }
   }

    /* n of type a^b - 1 */

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        mp_limb_t a, b, c, d, j, val, max_pow, base;
        mpz_t e, f, g, h;

        mpz_init(e);
        mpz_init(f);
        mpz_init(g);
        mpz_init(h);

        base = n_randint(state, upper_limit - 2) + 2;     /* base between 2 to 2642245*/
        max_pow = n_flog(UWORD_MAX, base);
        d = n_randint(state, max_pow);
        if (d < 2)                                /* root between 2 to max_pow */
            d = 2;

        flint_mpz_set_ui(h, d);
        c = n_pow(base, d) - 1;                   /* number */
        flint_mpz_set_ui(g, c);
        a = n_rootrem(&b, c, d);
        mpz_rootrem(e, f, g, flint_mpz_get_ui(h));

        val = flint_mpz_get_ui(e);
        j = flint_mpz_get_ui(f);

        result = ((a == val) && (b == j));

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("Passed Parameters : n = %wu root = %wu", c, d);
            flint_printf("Answer generated : base = %wu remainder = %wu", a, b);
            flint_printf("Expected answer : base = %wu remainder = %wu", val, j);
            fflush(stdout);
            flint_abort();
        }

        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
        mpz_clear(h);
    }

    TEST_FUNCTION_END(state);
}
