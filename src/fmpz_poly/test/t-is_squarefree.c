/*
    Copyright (C) 2009 William Hart
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

TEST_FUNCTION_START(fmpz_poly_is_squarefree, state)
{
    int i, result;

    /* Check that polynomials of degree <= 1 are square-free */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f;

        fmpz_poly_init(f);
        fmpz_poly_randtest(f, state, n_randint(state, 2), 100);

        result = (fmpz_poly_is_squarefree(f));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_debug(f), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
    }

    /* Check that a^2 f is not square-free */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, f;

        fmpz_poly_init(a);
        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1, 40);
        if (a->length < 2)
        {
            fmpz_poly_clear(a);
            continue;
        }
        fmpz_poly_init(f);
        fmpz_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1, 40);

        fmpz_poly_mul(a, a, a);
        fmpz_poly_mul(f, a, f);

        result = (!fmpz_poly_is_squarefree(f));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_debug(f), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(f);
    }

    /* Check that f + N*(x^M + 1) is square-free, for N >> f, M > deg(f) */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, f;
        fmpz_t N;

        fmpz_poly_init(a);
        fmpz_poly_set_coeff_si(a, 0, WORD(1));
        fmpz_poly_set_coeff_si(a, n_randint(state, 20), WORD(1));
        if (a->length < 2)
        {
            fmpz_poly_clear(a);
            continue;
        }
        fmpz_poly_init(f);
        fmpz_poly_randtest(f, state, a->length - 2, 40);

        fmpz_init_set_ui(N, UWORD(1));
        fmpz_mul_2exp(N, N, 45 + a->length);

        fmpz_poly_scalar_mul_fmpz(a, a, N);
        fmpz_poly_add(f, a, f);

        result = fmpz_poly_is_squarefree(f);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_debug(f), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(f);
        fmpz_clear(N);
    }

    TEST_FUNCTION_END(state);
}
