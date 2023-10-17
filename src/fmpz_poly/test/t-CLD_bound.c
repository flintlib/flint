/*
    Copyright (C) 2010, 2016 William Hart
    Copyright (C) 2010, Andy Novocin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_CLD_bound, state)
{
    int i, result = 1;

    /*
       test that CLD_bound is between the absolute value of the n-th
       coeff of f' and the sum of the absolute values of the coeffs of f'
    */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, g;
        fmpz_t sum, bound;
        slong xpow, j;
        slong bits = n_randint(state, 200) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);

        fmpz_init(sum);
        fmpz_init(bound);

        do
        {
           fmpz_poly_randtest(a, state, n_randint(state, 200), bits);

           if (!fmpz_poly_is_zero(a))
           {
              xpow = 0;

              while (fmpz_is_zero(a->coeffs + xpow))
                 xpow++;

              fmpz_poly_shift_right(a, a, xpow);
           }

           fmpz_poly_derivative(b, a);
           fmpz_poly_gcd(g, a, b);
        } while (!fmpz_poly_is_one(g));

        fmpz_poly_scalar_abs(b, b);

        for (j = 0; j < b->length && result; j++)
        {
           fmpz_poly_CLD_bound(bound, a, j);

           result &= (fmpz_cmp(b->coeffs + j, bound) <= 0);

           fmpz_add(sum, b->coeffs + j, b->coeffs + j);
        }

        result &= (fmpz_cmp(bound, sum) <= 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("length = %wd, bits = %wd, j = %wd\n", a->length, bits, j - 1);
            fmpz_print(b->coeffs + j - 1), flint_printf("\n\n");
            fmpz_print(bound), flint_printf("\n\n");
            fmpz_print(sum), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);

        fmpz_clear(sum);
        fmpz_clear(bound);
    }

    /*
       let f have a factor g (by setting f = g*h) then check that
       the bounds for f*g'/g are valid
    */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d, g;
        fmpz_t bound;
        slong xpow, j;
        slong bits = n_randint(state, 200) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_init(g);

        fmpz_init(bound);

        do
        {
           do {
              fmpz_poly_randtest(a, state, n_randint(state, 100) + 2, bits);

              if (!fmpz_poly_is_zero(a))
              {
                 xpow = 0;

                 while (fmpz_is_zero(a->coeffs + xpow))
                    xpow++;

                 fmpz_poly_shift_right(a, a, xpow);
              }
           } while (a->length < 2);

           do {
              fmpz_poly_randtest(b, state, n_randint(state, 100) + 1, bits);

              if (!fmpz_poly_is_zero(b))
              {
                 xpow = 0;

                 while (fmpz_is_zero(b->coeffs + xpow))
                    xpow++;

                 fmpz_poly_shift_right(b, b, xpow);
              }
           } while (fmpz_poly_is_zero(b));

           fmpz_poly_mul(c, a, b);

           fmpz_poly_derivative(d, c);
           fmpz_poly_gcd(g, c, d);
        } while (!fmpz_poly_is_one(g));

        fmpz_poly_derivative(g, a);
        fmpz_poly_mul(b, b, g);

        for (j = 0; j < b->length && result; j++)
        {
           fmpz_poly_CLD_bound(bound, c, j);

           result &= (fmpz_cmp(b->coeffs + j, bound) <= 0);
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("length = %wd, bits = %wd, j = %wd\n", a->length, bits, j - 1);
            fmpz_print(b->coeffs + j - 1), flint_printf("\n\n");
            fmpz_print(bound), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
        fmpz_poly_clear(g);

        fmpz_clear(bound);
    }

    TEST_FUNCTION_END(state);
}
