/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_mulhigh, state)
{
    int i, result;

    /* Compare with left truncated product of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        slong j, n;

        mp_limb_t m = n_randtest_not_zero(state);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        n = n_randint(state, 50);
        nmod_poly_randtest(b, state, n);
        nmod_poly_randtest(c, state, n);

        nmod_poly_mulhigh(a, b, c, n);
        nmod_poly_mul(b, b, c);
        for (j = 0; j < n; j++)
        {
            if (j < a->length)
                a->coeffs[j] = 0;
            if (j < b->length)
                b->coeffs[j] = 0;
        }
        _nmod_poly_normalise(a);
        _nmod_poly_normalise(b);

        result = (nmod_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
