/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly_q.h"

TEST_FUNCTION_START(fmpz_poly_q_set_equal, state)
{
    int i, result;

    /* Equal polynomials */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_q_t a, b;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_set(b, a);

        result = fmpz_poly_q_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_q_print(a), flint_printf("\n\n");
            fmpz_poly_q_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_q_t a, b;
        slong coeff = n_randint(state, 50);
        fmpz_t x1, x2;

        fmpz_init(x1);
        fmpz_init(x2);
        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_set(b, a);

        fmpz_poly_get_coeff_fmpz(x2, fmpz_poly_q_numref(b), coeff);
        do
            fmpz_randtest(x1, state, 50);
        while (fmpz_equal(x1, x2));
        fmpz_poly_set_coeff_fmpz(fmpz_poly_q_numref(b), coeff, x1);
        fmpz_poly_q_canonicalise(b);

        result = !fmpz_poly_q_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_q_print(a), flint_printf("\n\n");
            fmpz_poly_q_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x1);
        fmpz_clear(x2);
        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
    }

    TEST_FUNCTION_END(state);
}
