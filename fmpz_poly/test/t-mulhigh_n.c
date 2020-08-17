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
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("mulhigh_n....");
    fflush(stdout);

    

    /* Compare with left truncated product of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        slong j, n;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        n = n_randint(state, 50);
        fmpz_poly_randtest(b, state, n, 200);
        fmpz_poly_randtest(c, state, n, 200);

        fmpz_poly_mulhigh_n(a, b, c, n);
        fmpz_poly_mul(b, b, c);
        for (j = 0; j + 1 < n; j++)
        {
            if (j < a->length)
                fmpz_zero(a->coeffs + j);
            if (j < b->length)
                fmpz_zero(b->coeffs + j);
        }
        _fmpz_poly_normalise(a);
        _fmpz_poly_normalise(b);

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
        fmpz_poly_clear(c);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
