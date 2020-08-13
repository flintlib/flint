/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mod_poly.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("find_distinct_nonzero_roots....");
    fflush(stdout);

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        int highdegreefactor;
        fmpz_mod_poly_t a, b, r;
        fmpz_t p, e, zero;

        fmpz_init(zero);
        fmpz_init_set_ui(p, n_randtest_prime(state, 1));
        fmpz_init(e);
        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(r, p);

        fmpz_mod_poly_set_ui(a, 1 + n_randint(state, fmpz_get_ui(p) - 1));
        highdegreefactor = 0;
        for (j = n_randint(state, 10); j >= 0; j--)
        {
            if (n_randint(state, 10) > 1)
            {
                fmpz_mod_poly_randtest_monic_irreducible(b, state, 1);
            }
            else
            {
                highdegreefactor = 1;
                fmpz_mod_poly_randtest_monic_irreducible(b, state, 2 + n_randint(state, 9));
            }
            fmpz_mod_poly_mul(a, a, b);
        }

        fmpz_mod_poly_fit_length(r, fmpz_mod_poly_degree(a));
        if (fmpz_mod_poly_find_distinct_nonzero_roots(r->coeffs, a))
        {
            /* check that a is square free */
            fmpz_mod_poly_derivative(b, a);
            fmpz_mod_poly_gcd(b, b, a);
            if (fmpz_mod_poly_degree(b) > 0)
            {
                flint_printf("FAIL\ncheck multiple roots i = %wd\n", i);
                flint_abort();
            }

            /* check that each root is a root */
            for (j = fmpz_mod_poly_degree(a) - 1; j >= 0; j--)
            {
                if (fmpz_is_zero(r->coeffs + j))
                {
                    flint_printf("FAIL\ncheck zero root i = %wd\n", i);
                    flint_abort();
                }
                if (fmpz_mod_poly_evaluate_fmpz(e, a, r->coeffs + j), !fmpz_is_zero(e))
                {
                    flint_printf("FAIL\ncheck root is a root i = %wd\n", i);
                    flint_abort();
                }
            }
        }
        else
        {
            fmpz_mod_poly_derivative(b, a);
            fmpz_mod_poly_gcd(b, b, a);
            if (!highdegreefactor
                && fmpz_mod_poly_degree(b) == 0
                && (fmpz_mod_poly_evaluate_fmpz(e, a, zero), !fmpz_is_zero(e)))
            {
                flint_printf("FAIL\ncheck fail return i = %wd\n", i);
                flint_abort();
            }
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(r);
        fmpz_clear(p);
        fmpz_clear(e);
        fmpz_clear(zero);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
