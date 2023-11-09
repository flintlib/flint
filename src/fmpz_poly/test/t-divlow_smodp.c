/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_divlow_smodp, state)
{
    int i, result;

    /* Check a*b/a has correct coefficients */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong j, n;
        ulong p;
        fmpz_poly_t a, b, c, d;
        fmpz_t t, P, d1, d2;
        fmpz * vec;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_init(P);
        fmpz_init(t);
        fmpz_init(d1);
        fmpz_init(d2);

        p = n_randprime(state, n_randint(state, 7) + 2, 0);
        fmpz_set_ui(P, p);
        fmpz_pow_ui(P, P, n_randint(state, 4) + 1);

        n = n_randint(state, 40) + 2;

        vec = _fmpz_vec_init(n);

        do {
           fmpz_poly_randtest(a, state, n_randint(state, 100) + n, 200);
           if (a->length > 1)
           {
              fmpz_gcd(d1, P, a->coeffs + 0);
              fmpz_gcd(d2, P, a->coeffs + a->length - 1);
           }
        } while (a->length <= n || !fmpz_is_one(d1) || !fmpz_is_one(d2));

        do {
           fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
           if (b->length > 1)
           {
              fmpz_gcd(d1, P, b->coeffs + 0);
              fmpz_gcd(d2, P, b->coeffs + b->length - 1);
           }
        } while (b->length < 2 || !fmpz_is_one(d1) || !fmpz_is_one(d2));

        fmpz_poly_mul(c, a, b);

        fmpz_poly_scalar_mod_fmpz(d, b, P);

        fmpz_poly_divlow_smodp(vec, c, d, P, n);

        result = 1;
        for (j = 0; j < n && result; j++)
        {
           fmpz_smod(t, a->coeffs + j, P);
           result = fmpz_equal(t, vec + j);
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd, j = %wd\n\n", n, j);
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            fmpz_poly_print(d), flint_printf("\n\n");
            fmpz_print(P); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(d1);
        fmpz_clear(d2);
        fmpz_clear(t);
        fmpz_clear(P);
        _fmpz_vec_clear(vec, n);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
