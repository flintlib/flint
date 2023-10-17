/*
    Copyright (C) 2011 William Hart
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

TEST_FUNCTION_START(nmod_poly_invsqrt_series, state)
{
    int i, result;

    /* Check 1/g^2 = h mod x^m */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t h, g, r;
        slong m;

        mp_limb_t n;
        do n = n_randtest_prime(state, 0);
        while (n == UWORD(2));

        nmod_poly_init(h, n);
        nmod_poly_init(g, n);
        nmod_poly_init(r, n);

        do nmod_poly_randtest(h, state, n_randint(state, 1000));
        while (h->length == 0);
        nmod_poly_set_coeff_ui(h, 0, UWORD(1));

        m = n_randint(state, h->length) + 1;

        nmod_poly_invsqrt_series(g, h, m);

        nmod_poly_mullow(r, g, g, m);
        nmod_poly_inv_series(r, r, m);
        nmod_poly_truncate(h, m);

        result = (nmod_poly_equal(r, h));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(h), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(h);
        nmod_poly_clear(g);
        nmod_poly_clear(r);
    }

    /* Check aliasing of h and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t g, h;
        slong m;

        mp_limb_t n;
        do n = n_randtest_prime(state, 0);
        while (n == UWORD(2));

        nmod_poly_init(h, n);
        nmod_poly_init(g, n);
        do nmod_poly_randtest(h, state, n_randint(state, 500));
        while (h->length == 0);
        nmod_poly_set_coeff_ui(h, 0, UWORD(1));

        m = n_randint(state, h->length) + 1;

        nmod_poly_invsqrt_series(g, h, m);
        nmod_poly_invsqrt_series(h, h, m);

        result = (nmod_poly_equal(g, h));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(h), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("n = %wd, m = %wd\n", n, m);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(g);
        nmod_poly_clear(h);
    }

    TEST_FUNCTION_END(state);
}
