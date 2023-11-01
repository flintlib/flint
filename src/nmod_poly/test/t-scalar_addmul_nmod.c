/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_scalar_addmul_nmod, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t x = n_randint(state, n);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);

        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));
        nmod_poly_randtest(c, state, n_randint(state, 100));
        nmod_poly_randtest(d, state, n_randint(state, 100));

        nmod_poly_scalar_mul_nmod(c, b, x);
        nmod_poly_add(d, a, c);
        nmod_poly_scalar_addmul_nmod(a, b, x);
        if (!nmod_poly_equal(a, d))
        {
            flint_printf("FAIL: check non-aliasing\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));
        nmod_poly_randtest(c, state, n_randint(state, 100));
        nmod_poly_randtest(d, state, n_randint(state, 100));

        nmod_poly_scalar_mul_nmod(c, a, x);
        nmod_poly_add(d, a, c);
        nmod_poly_scalar_addmul_nmod(a, a, x);
        if (!nmod_poly_equal(a, d))
        {
            flint_printf("FAIL: check aliasing\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
