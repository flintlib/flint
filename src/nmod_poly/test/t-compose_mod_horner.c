/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_compose_mod_horner, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d, e;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(d, m);
        nmod_poly_init(e, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_compose_mod_horner(d, a, b, c);
        nmod_poly_compose(e, a, b);
        nmod_poly_rem(e, e, c);

        if (!nmod_poly_equal(d, e))
        {
            flint_printf("FAIL (composition):\n");
            nmod_poly_print(a); flint_printf("\n");
            nmod_poly_print(b); flint_printf("\n");
            nmod_poly_print(c); flint_printf("\n");
            nmod_poly_print(d); flint_printf("\n");
            nmod_poly_print(e); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
        nmod_poly_clear(e);
    }

    /* Test aliasing of res and a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(d, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_compose_mod_horner(d, a, b, c);
        nmod_poly_compose_mod_horner(a, a, b, c);

        if (!nmod_poly_equal(d, a))
        {
            flint_printf("FAIL (aliasing a):\n");
            nmod_poly_print(a); flint_printf("\n");
            nmod_poly_print(b); flint_printf("\n");
            nmod_poly_print(c); flint_printf("\n");
            nmod_poly_print(d); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    /* Test aliasing of res and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(d, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_compose_mod_horner(d, a, b, c);
        nmod_poly_compose_mod_horner(b, a, b, c);

        if (!nmod_poly_equal(d, b))
        {
            flint_printf("FAIL (aliasing b)\n");
            nmod_poly_print(a); flint_printf("\n");
            nmod_poly_print(b); flint_printf("\n");
            nmod_poly_print(c); flint_printf("\n");
            nmod_poly_print(d); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    /* Test aliasing of res and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t m = n_randtest_prime(state, 0);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(d, m);

        nmod_poly_randtest(a, state, 1+n_randint(state, 20));
        nmod_poly_randtest(b, state, 1+n_randint(state, 20));
        nmod_poly_randtest_not_zero(c, state, 1+n_randint(state, 20));

        nmod_poly_compose_mod_horner(d, a, b, c);
        nmod_poly_compose_mod_horner(c, a, b, c);

        if (!nmod_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing c)\n");
            nmod_poly_print(a); flint_printf("\n");
            nmod_poly_print(b); flint_printf("\n");
            nmod_poly_print(c); flint_printf("\n");
            nmod_poly_print(d); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
