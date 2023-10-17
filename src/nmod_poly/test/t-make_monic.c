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

TEST_FUNCTION_START(nmod_poly_make_monic, state)
{
    int i, result;

    /* Check new leading coeff = gcd old leading coeff and modulus */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t l;

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);

        if (n == 1) continue;
        do
        {
           nmod_poly_randtest(a, state, n_randint(state, 100) + 1);
        } while (a->length == 0 || n_gcd(*nmod_poly_lead(a), n) != 1);

        nmod_poly_make_monic(b, a);
        l = n_gcd(a->mod.n, a->coeffs[a->length - 1]);

        result = (l == b->coeffs[b->length - 1]);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("l = %wu, a->lead = %wd, n = %wu\n",
                l, a->coeffs[a->length - 1], a->mod.n);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    /* test aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t l;

        nmod_poly_init(a, n);

        if (n == 1) continue;
        do
        {
            nmod_poly_randtest(a, state, n_randint(state, 100) + 1);
        } while (a->length == 0 || n_gcd(*nmod_poly_lead(a), n) != 1);

        l = n_gcd(a->mod.n, a->coeffs[a->length - 1]);
        nmod_poly_make_monic(a, a);

        result = (l == a->coeffs[a->length - 1]);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("l = %wu, a->lead = %wd, n = %wu\n",
                l, a->coeffs[a->length - 1], a->mod.n);
            nmod_poly_print(a), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
    }

    TEST_FUNCTION_END(state);
}
