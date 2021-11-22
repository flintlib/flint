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

    flint_printf("xgcd....");
    fflush(stdout);

    

    /* Test aliasing of d and f, a and g */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, c, f, g, F, G;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(F);
        fmpz_init(G);

        fmpz_randtest_unsigned(G, state, 200);
        fmpz_add_ui(G, G, 1);
        fmpz_randm(F, state, G);
        if (n_randint(state, 2)) fmpz_neg(G, G);
        if (n_randint(state, 2)) fmpz_neg(F, F);
        fmpz_set(f, F);
        fmpz_set(g, G);

        fmpz_xgcd(d, a, b, f, g);
        fmpz_xgcd(f, g, c, f, g);

        result = (fmpz_equal(d, f)
               && fmpz_equal(b, c)
               && fmpz_equal(a, g));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("F = "), fmpz_print(F), flint_printf("\n");
            flint_printf("G = "), fmpz_print(G), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(F);
        fmpz_clear(G);
    }

    /* Test aliasing of a and f, d and g */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, c, f, g, F, G;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(F);
        fmpz_init(G);

        fmpz_randtest_unsigned(G, state, 200);
        fmpz_add_ui(G, G, 1);
        fmpz_randm(F, state, G);
        if (n_randint(state, 2)) fmpz_neg(G, G);
        if (n_randint(state, 2)) fmpz_neg(F, F);
        fmpz_set(f, F);
        fmpz_set(g, G);

        fmpz_xgcd(d, a, b, f, g);
        fmpz_xgcd(g, f, c, f, g);

        result = (fmpz_equal(d, g)
               && fmpz_equal(b, c)
               && fmpz_equal(a, f));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("F = "), fmpz_print(F), flint_printf("\n");
            flint_printf("G = "), fmpz_print(G), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(F);
        fmpz_clear(G);
    }

    /* Test aliasing of d and f, b and g */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, c, f, g, F, G;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(F);
        fmpz_init(G);

        fmpz_randtest_unsigned(G, state, 200);
        fmpz_add_ui(G, G, 1);
        fmpz_randm(F, state, G);
        if (n_randint(state, 2)) fmpz_neg(G, G);
        if (n_randint(state, 2)) fmpz_neg(F, F);
        fmpz_set(f, F);
        fmpz_set(g, G);

        fmpz_xgcd(d, a, b, f, g);
        fmpz_xgcd(f, c, g, f, g);

        result = (fmpz_equal(d, f) 
               && fmpz_equal(a, c)
               && fmpz_equal(b, g));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("F = "), fmpz_print(F), flint_printf("\n");
            flint_printf("G = "), fmpz_print(G), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(F);
        fmpz_clear(G);
    }

    /* Test aliasing of b and f, d and g */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, c, f, g, F, G;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(F);
        fmpz_init(G);

        fmpz_randtest_unsigned(G, state, 200);
        fmpz_add_ui(G, G, 1);
        fmpz_randm(F, state, G);
        if (n_randint(state, 2)) fmpz_neg(G, G);
        if (n_randint(state, 2)) fmpz_neg(F, F);
        fmpz_set(f, F);
        fmpz_set(g, G);

        fmpz_xgcd(d, a, b, f, g);
        fmpz_xgcd(g, c, f, f, g);

        result = (fmpz_equal(d, g)
               && fmpz_equal(a, c)
               && fmpz_equal(b, f));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("F = "), fmpz_print(F), flint_printf("\n");
            flint_printf("G = "), fmpz_print(G), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(F);
        fmpz_clear(G);
    }

    /* Test a f  + b g == d and d >= 0 */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, f, g, t1, t2;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(t1);
        fmpz_init(t2);

        fmpz_randtest_unsigned(g, state, 200);
        fmpz_add_ui(g, g, 1);
        fmpz_randm(f, state, g);
        if (n_randint(state, 2)) fmpz_neg(g, g);
        if (n_randint(state, 2)) fmpz_neg(f, f);
        
        fmpz_xgcd(d, a, b, f, g);

        fmpz_mul(t1, a, f);
        fmpz_mul(t2, b, g);
        fmpz_add(t1, t1, t2);

        result = fmpz_equal(t1, d) && fmpz_sgn(d) >= 0;
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("t1 = "), fmpz_print(t1), flint_printf("\n");
            flint_printf("t2 = "), fmpz_print(t2), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(t1);
        fmpz_clear(t2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

