/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <math.h>
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_get_mpf, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t x, y;
        mpz_t z;
        mpf_t a, b, tmp;

        fmpz_init(x);
        fmpz_init(y);
        mpz_init(z);
        mpf_inits(a, b, tmp, NULL);

        fmpz_randtest(x, state, 200);
        fmpz_get_mpz(z, x);

        fmpz_get_mpf(a, x);
        mpf_set_z(b, z);

        result = (mpf_cmp(a, b) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("a = "), mpf_out_str(stdout, 10, 0, a),
                flint_printf("\n");
            flint_printf("b = "), mpf_out_str(stdout, 10, 0, b),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        mpf_set_d(tmp,
                  (n_randtest(state) / (double) n_randtest_not_zero(state)));
        mpf_mul(a, a, tmp);

        fmpz_set_mpf(x, a);
        mpz_set_f(z, a);

        fmpz_set_mpz(y, z);
        result = fmpz_equal(x, y);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("y = "), fmpz_print(y), flint_printf("\n");
            flint_printf("a = "), mpf_out_str(stdout, 10, 0, a),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
        fmpz_clear(y);
        mpz_clear(z);
        mpf_clears(a, b, tmp, NULL);
    }

    TEST_FUNCTION_END(state);
}
