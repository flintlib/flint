/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("div_series....");
    fflush(stdout);

    /* Check A*B^{-1} * B is congruent A mod t^n */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, d;
        slong n = n_randint(state, 80) + 1;

        fmpz_init(p);
        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(d, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 80));
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 80) + 1);
        if (fmpz_is_zero(b->coeffs + 0))
           fmpz_add_ui(b->coeffs + 0, b->coeffs + 0, 1);
       
        fmpz_mod_poly_div_series(c, a, b, n);
        fmpz_mod_poly_mullow(d, c, b, n);

        result = (fmpz_mod_poly_equal_trunc(d, a, n));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("c = "), fmpz_mod_poly_print(c), flint_printf("\n\n");
            flint_printf("d = "), fmpz_mod_poly_print(d), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(d);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

