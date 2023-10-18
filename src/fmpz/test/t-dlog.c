/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_dlog, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t x;
        mpz_t z;
        mpfr_t r;
        double y, w;

        fmpz_init(x);
        mpz_init(z);
        mpfr_init2(r, 53);

        fmpz_randtest_not_zero(x, state, 10000);
        fmpz_abs(x, x);

        y = fmpz_dlog(x);

        fmpz_get_mpz(z, x);
        mpfr_set_z(r, z, MPFR_RNDN);

        mpfr_log(r, r, MPFR_RNDN);
        w = mpfr_get_d(r, MPFR_RNDN);

        result = (FLINT_ABS(y - w) <= w * 1e-13);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("y = %.20g\n", y);
            flint_printf("w = %.20g\n", w);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
        mpz_clear(z);
        mpfr_clear(r);
    }

    mpfr_free_cache();
    TEST_FUNCTION_END(state);
}
