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

TEST_FUNCTION_START(fmpz_poly_evaluate_mod, state)
{
    int i, result;

    /* Compare with evaluation over the integers */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t b, s;
        fmpz_poly_t f;
        mp_limb_t a, n, r;

        fmpz_poly_init(f);
        fmpz_poly_randtest(f, state, n_randint(state, 10), 20);

        n = n_randtest_not_zero(state);
        a = n_randint(state, n);

        fmpz_init(b);
        fmpz_init(s);
        fmpz_set_ui(b, a);

        r = fmpz_poly_evaluate_mod(f, a, n);
        fmpz_poly_evaluate_fmpz(s, f, b);

        result = (r == fmpz_mod_ui(s, s, n));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            gmp_printf("a = %Mu\n\n", a);
            gmp_printf("n = %Mu\n\n", n);
            gmp_printf("r = %Mu\n\n", r);
            flint_printf("s = "), fmpz_print(s), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_clear(b);
        fmpz_clear(s);
    }

    TEST_FUNCTION_END(state);
}
