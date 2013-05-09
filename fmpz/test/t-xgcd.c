/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009 William Hart

******************************************************************************/

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
    flint_rand_t state;

    printf("xgcd....");
    fflush(stdout);

    flint_randinit(state);

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
            printf("FAIL:\n\n");
            printf("d = "), fmpz_print(d), printf("\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            printf("f = "), fmpz_print(f), printf("\n");
            printf("g = "), fmpz_print(g), printf("\n");
            printf("F = "), fmpz_print(F), printf("\n");
            printf("G = "), fmpz_print(G), printf("\n");
            abort();
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
            printf("FAIL:\n\n");
            printf("d = "), fmpz_print(d), printf("\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            printf("f = "), fmpz_print(f), printf("\n");
            printf("g = "), fmpz_print(g), printf("\n");
            printf("F = "), fmpz_print(F), printf("\n");
            printf("G = "), fmpz_print(G), printf("\n");
            abort();
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
            printf("FAIL:\n\n");
            printf("d = "), fmpz_print(d), printf("\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            printf("f = "), fmpz_print(f), printf("\n");
            printf("g = "), fmpz_print(g), printf("\n");
            printf("F = "), fmpz_print(F), printf("\n");
            printf("G = "), fmpz_print(G), printf("\n");
            abort();
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
            printf("FAIL:\n\n");
            printf("d = "), fmpz_print(d), printf("\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            printf("f = "), fmpz_print(f), printf("\n");
            printf("g = "), fmpz_print(g), printf("\n");
            printf("F = "), fmpz_print(F), printf("\n");
            printf("G = "), fmpz_print(G), printf("\n");
            abort();
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

    /* Test a f  + b g == d */
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

        result = fmpz_equal(t1, d);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("d = "), fmpz_print(d), printf("\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("f = "), fmpz_print(f), printf("\n");
            printf("g = "), fmpz_print(g), printf("\n");
            printf("t1 = "), fmpz_print(t1), printf("\n");
            printf("t2 = "), fmpz_print(t2), printf("\n");
            abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(t1);
        fmpz_clear(t2);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

