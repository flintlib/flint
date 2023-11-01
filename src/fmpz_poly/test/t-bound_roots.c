/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_bound_roots, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t f, g;
        fmpz_t t, p, q, bound, nbound;
        slong i, num_roots;

        fmpz_init(t);
        fmpz_init(p);
        fmpz_init(q);
        fmpz_init(bound);
        fmpz_init(nbound);
        fmpz_poly_init(f);
        fmpz_poly_init(g);

        /* start with a constant */
        fmpz_poly_randtest(f, state, 1, 1 + n_randint(state, 200));
        fmpz_zero(bound);

        num_roots = n_randint(state, 10);
        if (fmpz_poly_is_zero(f))
            num_roots = 0;

        for (i = 0; i < num_roots; i++)
        {
            fmpz_randtest(p, state, 1 + n_randint(state, 200));
            fmpz_randtest_not_zero(q, state, 1 + n_randint(state, 200));

            fmpz_abs(p, p);
            fmpz_abs(q, q);

            fmpz_cdiv_q(t, p, q);
            if (fmpz_cmp(t, bound) > 0)
                fmpz_set(bound, t);

            if (n_randint(state, 2))
                fmpz_neg(p, p);

            fmpz_poly_set_coeff_fmpz(g, 0, p);
            fmpz_poly_set_coeff_fmpz(g, 1, q);

            fmpz_poly_mul(f, f, g);
        }

        fmpz_poly_bound_roots(nbound, f);

        if (fmpz_cmp(nbound, bound) < 0)
        {
            flint_printf("FAIL\n");
            flint_printf("f = "); fmpz_poly_print(f); flint_printf("\n\n");
            flint_printf("bound = "); fmpz_print(bound); flint_printf("\n\n");
            flint_printf("computed bound = "); fmpz_print(nbound); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
       }

        fmpz_clear(t);
        fmpz_clear(p);
        fmpz_clear(q);
        fmpz_clear(bound);
        fmpz_clear(nbound);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
