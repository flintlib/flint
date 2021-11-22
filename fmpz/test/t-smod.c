/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Abhinav Baid

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

    flint_printf("smod....");
    fflush(stdout);

    for (i = 0; i < 200000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d, e;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(e);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        if (n_randint(state, 10) == 0)
        {
            fmpz_fdiv_q_ui(a, b, 2);
            if (n_randint(state, 2))
                fmpz_add(a, b, b);
            if (n_randint(state, 2))
                fmpz_sub(a, b, b);
            if (n_randint(state, 2))
                fmpz_add(a, b, b);
        }

        fmpz_smod(c, a, b);

        fmpz_sub(d, a, c);
        fmpz_mod(d, d, b);
        if (!fmpz_is_zero(d))
        {
            flint_printf("FAIL: check b|(smod(a,b) - a)\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("smod(a,b) = "), fmpz_print(c), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_abs(e, b);
        fmpz_mul_2exp(d, c, 1);
        if (fmpz_cmp(d, e) > 0)
        {
            flint_printf("FAIL: check 2*smod(a,b) <= |b|\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("smod(a,b) = "), fmpz_print(c), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_neg(e, e);
        if (fmpz_cmp(e, d) >= 0)
        {
            flint_printf("FAIL: check -|b| < 2*smod(a,b)\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("smod(a,b) = "), fmpz_print(c), flint_printf("\n");
            fflush(stdout);
            flint_abort();
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

        fmpz_smod(c, a, a);
        fmpz_smod(d, a, b);

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
            fflush(stdout);
            flint_abort();
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

        fmpz_smod(a, a, b);
        fmpz_smod(d, c, b);

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
            fflush(stdout);
            flint_abort();
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

        fmpz_smod(b, a, b);
        fmpz_smod(d, a, c);

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
            fflush(stdout);
            flint_abort();
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
