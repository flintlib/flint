/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("get_mpfr....");
    fflush(stdout);

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
            abort();
        }

        fmpz_clear(x);
        fmpz_clear(y);
        mpz_clear(z);
        mpfr_clears(a, b, NULL);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
