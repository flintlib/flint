/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("scalar_mul_si....");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        slong n = z_randtest(state);

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpz_poly_scalar_mul_si(b, a, n);
        fmpz_poly_scalar_mul_si(a, a, n);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* Compare with fmpz_poly_scalar_mul_ui */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        ulong n = n_randbits(state, FLINT_BITS - 1);

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpz_poly_scalar_mul_ui(b, a, n);
        fmpz_poly_scalar_mul_si(a, a, n);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* Check (a*n1)*n2 = a*(n1*n2) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        slong n1 = (slong) n_randbits(state, (FLINT_BITS - 2) / 2);
        slong n2 = (slong) n_randbits(state, (FLINT_BITS - 2) / 2);
        if (n_randint(state, 2))
            n1 = -n1;
        if (n_randint(state, 2))
            n2 = -n2;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpz_poly_scalar_mul_si(b, a, n1);
        fmpz_poly_scalar_mul_si(c, b, n2);
        fmpz_poly_scalar_mul_si(b, a, n1 * n2);

        result = (fmpz_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL n1 = %wd, n2 = %wd:\n", n1, n2);
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
