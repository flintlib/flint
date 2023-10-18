/*
    Copyright (C) 2009, 2013 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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

TEST_FUNCTION_START(fmpz_poly_rem_powers_precomp, state)
{
    int i, result;

    /* Compare with full division, no aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q, r, r2;
        fmpz_poly_powers_precomp_t b_inv;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_init(r2);

        fmpz_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, a->length + 1) + 1, 100);
        fmpz_set_ui(b->coeffs + b->length - 1, 1); /* b must be monic */

        fmpz_poly_divrem_basecase(q, r, a, b);
        fmpz_poly_powers_precompute(b_inv, b);
        fmpz_poly_rem_powers_precomp(r2, a, b, b_inv);

        result = (fmpz_poly_equal(r, r2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            fmpz_poly_print(r2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_powers_clear(b_inv);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
        fmpz_poly_clear(r2);
    }

    /* Check q and a alias */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, b_inv, q;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(b_inv);
        fmpz_poly_init(q);
        fmpz_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);
        fmpz_set_ui(b->coeffs + b->length - 1, 1); /* b must be monic */

        fmpz_poly_div_basecase(q, a, b);
        fmpz_poly_preinvert(b_inv, b);
        fmpz_poly_div_preinv(a, a, b, b_inv);

        result = (fmpz_poly_equal(a, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(b_inv);
        fmpz_poly_clear(q);
    }

    /* Check q and b alias */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, b_inv, q;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(b_inv);
        fmpz_poly_init(q);
        fmpz_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);
        fmpz_set_ui(b->coeffs + b->length - 1, 1); /* b must be monic */

        fmpz_poly_div_basecase(q, a, b);
        fmpz_poly_preinvert(b_inv, b);
        fmpz_poly_div_preinv(b, a, b, b_inv);

        result = (fmpz_poly_equal(b, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(b_inv);
        fmpz_poly_clear(q);
    }

    TEST_FUNCTION_END(state);
}
