/*
    Copyright (C) 2010, 2020 William Hart

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

    flint_printf("mullow_SS_precache....");
    fflush(stdout);

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
	fmpz_poly_mul_precache_t pre;
        slong len, trunc;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

        len = b->length + c->length - 1;
        trunc = (len <= 0) ? 0 : n_randint(state, b->length + c->length);

        fmpz_poly_mul_SS_precache_init(pre, 50, 200, c);
	
        fmpz_poly_mullow_SS(a, b, c, trunc);
        fmpz_poly_mullow_SS_precache(b, b, pre, trunc);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            abort();
        }

	fmpz_poly_mul_precache_clear(pre);

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Compare with mul_KS */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d;
        fmpz_poly_mul_precache_t pre;
	slong len, trunc;
	int k;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);
	
	fmpz_poly_mul_SS_precache_init(pre, 50, 200, c);
	
	for (k = 0; k < 3; k++)
	{
	   fmpz_poly_randtest(b, state, n_randint(state, 50), 200);

           len = b->length + c->length - 1;
           trunc = (len <= 0) ? 0 : n_randint(state, b->length + c->length - 1);

           fmpz_poly_mul_KS(a, b, c);
           fmpz_poly_truncate(a, trunc);
           fmpz_poly_mullow_SS_precache(d, b, pre, trunc);

           result = (fmpz_poly_equal(a, d));
           if (!result)
           {
               flint_printf("FAIL:\n");
               fmpz_poly_print(a), flint_printf("\n\n");
               fmpz_poly_print(d), flint_printf("\n\n");
               abort();
           }
	}

	fmpz_poly_mul_precache_clear(pre);

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
