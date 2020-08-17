/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

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

    flint_printf("revert_series_lagrange_fast....");
    fflush(stdout);

    

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g;
        slong n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_randtest(g, state, n_randint(state, 50),
            1+n_randint(state,100));
        fmpz_poly_set_coeff_ui(g, 0, 0);
        fmpz_poly_set_coeff_ui(g, 1, 1);
        if (n_randlimb(state) % 2)
            fmpz_poly_neg(g, g);  /* get -x term */
        n = n_randint(state, 50);

        fmpz_poly_revert_series_lagrange_fast(f, g, n);
        fmpz_poly_revert_series_lagrange_fast(g, g, n);

        result = (fmpz_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            fmpz_poly_print(g), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    /* Check f(f^(-1)) = id */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h;
        slong n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 50), 1+n_randint(state,100));
        fmpz_poly_set_coeff_ui(g, 0, 0);
        fmpz_poly_set_coeff_ui(g, 1, 1);
        if (n_randlimb(state) % 2)
            fmpz_poly_neg(g, g);  /* get -x term */
        n = n_randint(state, 50);

        fmpz_poly_revert_series_lagrange_fast(f, g, n);
        fmpz_poly_compose_series(h, g, f, n);

        result = ((n <= 1 && fmpz_poly_is_zero(h)) ||
            (h->length == 2 && fmpz_is_zero(h->coeffs + 0) &&
                fmpz_is_one(h->coeffs + 1)));
        if (!result)
        {
            flint_printf("FAIL (comparison):\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            fmpz_poly_print(g), flint_printf("\n\n");
            fmpz_poly_print(h), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
