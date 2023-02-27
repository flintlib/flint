/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2009 Andy Novocin

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

    double d, d2;
    slong exp, exp2;

    FLINT_TEST_INIT(state);

    flint_printf("set_d_2exp....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        fmpz_init(a);

        fmpz_randtest(a, state, 200);

        d = fmpz_get_d_2exp(&exp, a);
        fmpz_set_d_2exp(a, d, exp);
        d2 = fmpz_get_d_2exp(&exp2, a);

        result = (d2 == d && exp == exp2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("d = %f, d2 = %f\n", d, d2);
            flint_printf("exp = %wd, exp2 = %wd, bits = %wu\n", exp, exp2, fmpz_bits(a));
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
