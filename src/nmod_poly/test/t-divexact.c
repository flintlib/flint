/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_divexact, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, ab, q;
        int aliasing;

        ulong n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(ab, n);
        nmod_poly_init(q, n);

        nmod_poly_randtest(a, state, n_randint(state, 400));
        do nmod_poly_randtest(b, state, n_randint(state, 400));
        while (b->length == 0);

        nmod_poly_mul(ab, a, b);

        aliasing = n_randint(state, 3);

        switch (aliasing)
        {
            case 0:
                nmod_poly_divexact(q, ab, b);
                break;
            case 1:
                nmod_poly_set(q, ab);
                nmod_poly_divexact(q, q, b);
                break;
            default:
                nmod_poly_set(q, b);
                nmod_poly_divexact(q, ab, q);
                break;
        }

        result = (nmod_poly_equal(q, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(ab);
        nmod_poly_clear(q);
    }

    TEST_FUNCTION_END(state);
}
