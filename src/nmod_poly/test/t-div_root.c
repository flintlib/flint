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

TEST_FUNCTION_START(nmod_poly_div_root, state)
{
    int i, result;

    /* Compare with standard divrem */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t P, Q, D, DQ, DR;
        mp_limb_t mod, r, rem;
        slong n;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 100);
        r = n_randint(state, mod);

        nmod_poly_init(P, mod);
        nmod_poly_init(Q, mod);
        nmod_poly_init(D, mod);
        nmod_poly_init(DQ, mod);
        nmod_poly_init(DR, mod);

        nmod_poly_randtest(P, state, n);

        rem = nmod_poly_div_root(Q, P, r);

        nmod_poly_set_coeff_ui(D, 0, n_negmod(r, mod));
        nmod_poly_set_coeff_ui(D, 1, UWORD(1));

        nmod_poly_divrem(DQ, DR, P, D);

        result = nmod_poly_equal(Q, DQ) &&
            (rem == nmod_poly_get_coeff_ui(DR, 0));

        if (!result)
        {
            flint_printf("FAIL!\n");
            flint_printf("P:\n"); nmod_poly_print(P); flint_printf("\n\n");
            flint_printf("Q:\n"); nmod_poly_print(Q); flint_printf("\n\n");
            flint_printf("D:\n"); nmod_poly_print(D); flint_printf("\n\n");
            flint_printf("DQ:\n"); nmod_poly_print(DQ); flint_printf("\n\n");
            flint_printf("DR:\n"); nmod_poly_print(DR); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(P);
        nmod_poly_clear(Q);
        nmod_poly_clear(D);
        nmod_poly_clear(DQ);
        nmod_poly_clear(DR);
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t P, Q1, Q2;
        mp_limb_t mod, r, rem1, rem2;
        slong n;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 100);
        r = n_randint(state, mod);

        nmod_poly_init(P, mod);
        nmod_poly_init(Q1, mod);
        nmod_poly_init(Q2, mod);

        nmod_poly_randtest(P, state, n);
        nmod_poly_set(Q2, P);

        rem1 = nmod_poly_div_root(Q1, P, r);
        rem2 = nmod_poly_div_root(Q2, Q2, r);

        result = nmod_poly_equal(Q1, Q2) && (rem1 == rem2);

        if (!result)
        {
            flint_printf("FAIL (aliasing)!\n");
            flint_printf("P:\n"); nmod_poly_print(P); flint_printf("\n\n");
            flint_printf("Q1:\n"); nmod_poly_print(Q1); flint_printf("\n\n");
            flint_printf("Q2:\n"); nmod_poly_print(Q2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(P);
        nmod_poly_clear(Q1);
        nmod_poly_clear(Q2);
    }

    TEST_FUNCTION_END(state);
}
