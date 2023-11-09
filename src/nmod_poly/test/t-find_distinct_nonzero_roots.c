/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_find_distinct_nonzero_roots, state)
{
    slong i, j;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        int highdegreefactor;
        nmod_poly_t a, b, r;
        mp_limb_t p;

        p = n_randtest_prime(state, 1);

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(r, p);

        nmod_poly_zero(a);
        nmod_poly_set_coeff_ui(a, 0, 1 + n_randint(state, p - 1));
        highdegreefactor = 0;
        for (j = n_randint(state, 10); j >= 0; j--)
        {
            if (n_randint(state, 10) > 1)
            {
                nmod_poly_randtest_monic_irreducible(b, state, 1);
            }
            else
            {
                highdegreefactor = 1;
                nmod_poly_randtest_monic_irreducible(b, state, 2 + n_randint(state, 9));
            }
            nmod_poly_mul(a, a, b);
        }

        nmod_poly_fit_length(r, nmod_poly_degree(a));
        if (nmod_poly_find_distinct_nonzero_roots(r->coeffs, a))
        {
            /* check that a is square free */
            nmod_poly_derivative(b, a);
            nmod_poly_gcd(b, b, a);
            if (nmod_poly_degree(b) > 0)
            {
                flint_printf("FAIL\ncheck multiple roots i = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }

            /* check that each root is a root */
            for (j = nmod_poly_degree(a) - 1; j >= 0; j--)
            {
                if (r->coeffs[j] == 0)
                {
                    flint_printf("FAIL\ncheck zero root i = %wd\n", i);
                    fflush(stdout);
                    flint_abort();
                }
                if (nmod_poly_evaluate_nmod(a, r->coeffs[j]) != 0)
                {
                    flint_printf("FAIL\ncheck root is a root i = %wd\n", i);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }
        else
        {
            nmod_poly_derivative(b, a);
            nmod_poly_gcd(b, b, a);
            if (!highdegreefactor
                && nmod_poly_degree(b) == 0
                && nmod_poly_evaluate_nmod(a, 0) != 0)
            {
                flint_printf("FAIL\ncheck fail return i = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(r);
    }

    TEST_FUNCTION_END(state);
}
