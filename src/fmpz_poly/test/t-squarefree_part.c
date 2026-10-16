/*
    Copyright (C) 2017 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_squarefree_part, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t p, p2, q1, q2;

        fmpz_poly_init(p);
        fmpz_poly_init(p2);
        fmpz_poly_init(q1);
        fmpz_poly_init(q2);

        fmpz_poly_randtest(p, state, 1 + n_randint(state, 20), 10 + n_randint(state, 100));

        fmpz_poly_squarefree_part(q1, p);

        if (fmpz_poly_is_squarefree(p))
        {
            fmpz_poly_primitive_part(q2, p);
            if (!fmpz_poly_equal(q1, q2))
            {
                flint_printf("FAIL:\n");
                flint_printf("q2 = "); fmpz_poly_print(q2); flint_printf("\n");
                flint_printf("q1 = "); fmpz_poly_print(q1); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_poly_mul(p2, p, p);
        fmpz_poly_squarefree_part(q2, p2);
        if (!fmpz_poly_equal(q1, q2))
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "); fmpz_poly_print(p); flint_printf("\n");
            flint_printf("q1 = "); fmpz_poly_print(q1); flint_printf("\n");
            flint_printf("q2 = "); fmpz_poly_print(q2); flint_printf("\n");
            flint_abort();
        }

        fmpz_poly_clear(p);
        fmpz_poly_clear(p2);
        fmpz_poly_clear(q1);
        fmpz_poly_clear(q2);
    }

    TEST_FUNCTION_END(state);
}
