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
#include <stdio.h>
#include <mpfr.h>
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_get_mpfr, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t x, y;
        mpz_t z;
        mpfr_t a, b;

        fmpz_init(x);
        fmpz_init(y);
        mpz_init(z);
        mpfr_inits(a, b, NULL);

        fmpz_randtest(x, state, 200);
        fmpz_get_mpz(z, x);

        fmpz_get_mpfr(a, x, MPFR_RNDN);
        mpfr_set_z(b, z, MPFR_RNDN);

        result = (mpfr_equal_p(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("a = "), mpfr_out_str(stdout, 10, 0, a, MPFR_RNDN),
                flint_printf("\n");
            flint_printf("b = "), mpfr_out_str(stdout, 10, 0, b, MPFR_RNDN),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
        fmpz_clear(y);
        mpz_clear(z);
        mpfr_clears(a, b, NULL);
    }

    TEST_FUNCTION_END(state);
}
