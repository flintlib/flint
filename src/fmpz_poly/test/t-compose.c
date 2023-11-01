/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2016 Shivin Srivastava

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_compose, state)
{
    int i, result;

    /* Bill's bug */
    {
        fmpz_poly_t f, g, h, s, t;
        slong k;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(s);
        fmpz_poly_init(t);

        fmpz_poly_set_str(g, "21  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6");
        fmpz_poly_set_str(h, "8  -2411740686 -274861162464 -4294966273 -35167058005888 4261511676 -1 8589869056 -70334401183747");

        fmpz_poly_set_ui(t, 1);
        for (k = 0; k < g->length; k++)
        {
            fmpz_poly_scalar_addmul_fmpz(s, t, g->coeffs + k);
            fmpz_poly_mul(t, t, h);
        }

        fmpz_poly_compose(f, g, h);

        result = (fmpz_poly_equal(f, s));
        if (!result)
        {
            flint_printf("FAIL (Bill's bug):\n");
            flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n\n");
            flint_printf("h = "), fmpz_poly_print(h), flint_printf("\n\n");
            flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
            flint_printf("s = "), fmpz_poly_print(s), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    /* Check aliasing of the first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpz_poly_compose(f, g, h);
        fmpz_poly_compose(g, g, h);

        result = (fmpz_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL (aliasing 1):\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            fmpz_poly_print(g), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpz_poly_compose(f, g, h);
        fmpz_poly_compose(h, g, h);

        result = (fmpz_poly_equal(f, h));
        if (!result)
        {
            flint_printf("FAIL (aliasing 2):\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            fmpz_poly_print(h), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    /* Compare with the naive method for g(h(t)) */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h, s, t;
        slong k;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(s);
        fmpz_poly_init(t);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpz_poly_set_ui(t, 1);
        for (k = 0; k < g->length; k++)
        {
            fmpz_poly_scalar_addmul_fmpz(s, t, g->coeffs + k);
            fmpz_poly_mul(t, t, h);
        }

        fmpz_poly_compose(f, g, h);

        result = (fmpz_poly_equal(f, s));
        if (!result)
        {
            flint_printf("FAIL (comparison):\n");
            flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n\n");
            flint_printf("h = "), fmpz_poly_print(h), flint_printf("\n\n");
            flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
            flint_printf("s = "), fmpz_poly_print(s), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    /* Testing linear composition*/
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h, s, t;
        fmpz_t c;
        slong k;

        fmpz_init(c);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(s);
        fmpz_poly_init(t);

        /* h is a linear polynomial*/
        for (k = 0; k < 2; k++)
        {
            fmpz_randtest(c, state, 50);
            fmpz_poly_set_coeff_fmpz(h, k, c);
        }

        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_set_ui(t, 1);
        for (k = 0; k < g->length; k++)
        {
            fmpz_poly_scalar_addmul_fmpz(s, t, g->coeffs + k);
            fmpz_poly_mul(t, t, h);
        }

        fmpz_poly_compose(f, g, h);

        result = (fmpz_poly_equal(f, s));
        if (!result)
        {
            flint_printf("FAIL (linear composition):\n");
            flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n\n");
            flint_printf("h = "), fmpz_poly_print(h), flint_printf("\n\n");
            flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
            flint_printf("s = "), fmpz_poly_print(s), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
