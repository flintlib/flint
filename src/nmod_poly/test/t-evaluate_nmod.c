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

TEST_FUNCTION_START(nmod_poly_evaluate_nmod, state)
{
    int i, j, result = 1;

    /* Check evaluation at 1 gives sum of coeffs */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t sum, eval;

        nmod_poly_init(a, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));

        eval = nmod_poly_evaluate_nmod(a, 1);

        sum = 0;
        for (j = 0; j < a->length; j++)
           sum = n_addmod(sum, nmod_poly_get_coeff_ui(a, j), n);

        result = (sum == eval);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu\n", a->length, a->mod.n);
            flint_printf("sum = %wu, eval = %wu\n", sum, eval);
            nmod_poly_print(a), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
    }

    /* Check a(c) + b(c) = (a + b)(c) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t eval1, eval2, c;

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        c = n_randint(state, n);

        eval1 = nmod_poly_evaluate_nmod(a, c);
        eval1 = n_addmod(eval1, nmod_poly_evaluate_nmod(b, c), n);

        nmod_poly_add(a, a, b);
        eval2 = nmod_poly_evaluate_nmod(a, c);

        result = (eval1 == eval2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("eval1 = %wu, eval2 = %wu\n", eval1, eval2);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
