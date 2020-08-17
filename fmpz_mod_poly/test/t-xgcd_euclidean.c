/*
    Copyright (C) 2012 Sebastian Pancratz

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
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("xgcd_euclidean....");
    fflush(stdout);

    

    /* Generic case, most likely co-prime arguments ******************************/

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, d, g, s, t, v, w;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_init(v, p);
        fmpz_mod_poly_init(w, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));

        fmpz_mod_poly_gcd_euclidean(d, a, b);
        fmpz_mod_poly_xgcd_euclidean(g, s, t, a, b);

        fmpz_mod_poly_mul(v, s, a);
        fmpz_mod_poly_mul(w, t, b);
        fmpz_mod_poly_add(w, v, w);

        result = (fmpz_mod_poly_equal(d, g) && fmpz_mod_poly_equal(g, w));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a), flint_printf("\n\n");
            fmpz_mod_poly_print(b), flint_printf("\n\n");
            fmpz_mod_poly_print(d), flint_printf("\n\n");
            fmpz_mod_poly_print(g), flint_printf("\n\n");
            fmpz_mod_poly_print(s), flint_printf("\n\n");
            fmpz_mod_poly_print(t), flint_printf("\n\n");
            fmpz_mod_poly_print(v), flint_printf("\n\n");
            fmpz_mod_poly_print(w), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(d);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(t);
        fmpz_mod_poly_clear(v);
        fmpz_mod_poly_clear(w);
        fmpz_clear(p);
    }

    /* Special case, arguments share a factor ********************************/

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, d, f, g, s, t, v, w;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(t, p);
        fmpz_mod_poly_init(v, p);
        fmpz_mod_poly_init(w, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(f, state, n_randint(state, 20));
        fmpz_mod_poly_mul(a, a, f);
        fmpz_mod_poly_mul(b, b, f);

        fmpz_mod_poly_gcd_euclidean(d, a, b);
        fmpz_mod_poly_xgcd_euclidean(g, s, t, a, b);

        fmpz_mod_poly_mul(v, s, a);
        fmpz_mod_poly_mul(w, t, b);
        fmpz_mod_poly_add(w, v, w);

        result = (fmpz_mod_poly_equal(d, g) && fmpz_mod_poly_equal(g, w));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a), flint_printf("\n\n");
            fmpz_mod_poly_print(b), flint_printf("\n\n");
            fmpz_mod_poly_print(d), flint_printf("\n\n");
            fmpz_mod_poly_print(f), flint_printf("\n\n");
            fmpz_mod_poly_print(g), flint_printf("\n\n");
            fmpz_mod_poly_print(s), flint_printf("\n\n");
            fmpz_mod_poly_print(t), flint_printf("\n\n");
            fmpz_mod_poly_print(v), flint_printf("\n\n");
            fmpz_mod_poly_print(w), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(d);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(t);
        fmpz_mod_poly_clear(v);
        fmpz_mod_poly_clear(w);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

