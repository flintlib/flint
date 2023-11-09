/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_pow_trunc_binexp, state)
{
    int i, result;

    /* Check powering against naive method */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        mp_limb_t n = n_randtest_not_zero(state);
        slong e, trunc;

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(a, state, n_randint(state, 30));
        e = n_randint(state, 20);
        trunc = n_randint(state, 30);

        nmod_poly_pow_trunc_binexp(b, a, e, trunc);

        nmod_poly_pow(c, a, e);
        nmod_poly_truncate(c, trunc);

        result = (nmod_poly_equal(b, c)
            || (a->length == 0 && e == 0 && c->length == 1 && c->coeffs[0] == 1));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu, exp = %wd, trunc = %wd\n",
                a->length, a->mod.n, e, trunc);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        mp_limb_t n = n_randtest_not_zero(state);
        slong e, trunc;

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(a, state, n_randint(state, 30));
        e = n_randint(state, 20);
        trunc = n_randint(state, 30);

        nmod_poly_pow_trunc_binexp(b, a, e, trunc);

        nmod_poly_set(c, a);
        nmod_poly_pow_trunc_binexp(c, c, e, trunc);

        result = (nmod_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu, exp = %wd, trunc = %wd\n",
                a->length, a->mod.n, e, trunc);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
