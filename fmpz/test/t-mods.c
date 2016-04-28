/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
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

    flint_printf("mods....");
    fflush(stdout);



    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d, e;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(e);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        fmpz_mods(c, a, b);
        fmpz_sub(d, a, c);
        fmpz_mod(d, d, b);
        fmpz_abs(e, b);
        fmpz_fdiv_q_2exp(e, e, 1);

        result = (fmpz_is_zero(d) && fmpz_cmp(c, e) <= 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = ");
            fmpz_print(a);
            flint_printf("\n");
            flint_printf("b = ");
            fmpz_print(b);
            flint_printf("\n");
            flint_printf("c = ");
            fmpz_print(c);
            flint_printf("\n");
            flint_printf("d = ");
            fmpz_print(d);
            flint_printf("\n");
            flint_printf("e = ");
            fmpz_print(e);
            flint_printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);

        fmpz_randtest_not_zero(a, state, 200);
        fmpz_set(b, a);

        fmpz_mods(c, a, a);
        fmpz_mods(d, a, b);

        result = (fmpz_cmp(c, d) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = ");
            fmpz_print(a);
            flint_printf("\n");
            flint_printf("b = ");
            fmpz_print(b);
            flint_printf("\n");
            flint_printf("c = ");
            fmpz_print(c);
            flint_printf("\n");
            flint_printf("d = ");
            fmpz_print(d);
            flint_printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
    }

    /* Test aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);
        fmpz_set(c, a);

        fmpz_mods(a, a, b);
        fmpz_mods(d, c, b);

        result = (fmpz_cmp(a, d) == 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = ");
            fmpz_print(a);
            flint_printf("\n");
            flint_printf("b = ");
            fmpz_print(b);
            flint_printf("\n");
            flint_printf("c = ");
            fmpz_print(c);
            flint_printf("\n");
            flint_printf("d = ");
            fmpz_print(d);
            flint_printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
    }

    /* Test aliasing of b and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);
        fmpz_set(c, b);

        fmpz_mods(b, a, b);
        fmpz_mods(d, a, c);

        result = (fmpz_cmp(b, d) == 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = ");
            fmpz_print(a);
            flint_printf("\n");
            flint_printf("b = ");
            fmpz_print(b);
            flint_printf("\n");
            flint_printf("c = ");
            fmpz_print(c);
            flint_printf("\n");
            flint_printf("d = ");
            fmpz_print(d);
            flint_printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
