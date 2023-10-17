/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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

/*
    Tests whether the polynomial is suitably normalised for the
    result of a GCD operation, that is, whether it's leading
    coefficient is non-negative.
 */
#ifndef _t_gcd_is_canonical
#define _t_gcd_is_canonical _t_gcd_is_canonical
static
int _t_gcd_is_canonical(const fmpz_poly_t poly)
{
    return fmpz_poly_is_zero(poly) || (fmpz_sgn(fmpz_poly_lead(poly)) > 0);
}
#endif

TEST_FUNCTION_START(fmpz_poly_gcd_heuristic, state)
{
    int i, result, d1, d2;

    /* Check aliasing of a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        d1 = fmpz_poly_gcd_heuristic(a, b, c);
        d2 = fmpz_poly_gcd_heuristic(b, b, c);

        result = ((d1 == 0 && d2 == 0) || (fmpz_poly_equal(a, b)
                                           && _t_gcd_is_canonical(a)));
        if (!result)
        {
            flint_printf("FAIL (aliasing a and b):\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        d1 = fmpz_poly_gcd_heuristic(a, b, c);
        d2 = fmpz_poly_gcd_heuristic(c, b, c);

        result = ((d1 == 0 && d2 == 0) || (fmpz_poly_equal(a, c)
                                           && _t_gcd_is_canonical(a)));
        if (!result)
        {
            flint_printf("FAIL (aliasing a and c):\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("c = "), fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check that a divides GCD(af, ag) */
    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, d, f, g, q, r;

        fmpz_poly_init(a);
        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 100) + 1, 40);
        fmpz_poly_randtest(f, state, n_randint(state, 100), 40);
        fmpz_poly_randtest(g, state, n_randint(state, 100), 40);

        fmpz_poly_mul(f, a, f);
        fmpz_poly_mul(g, a, g);
        d1 = fmpz_poly_gcd_heuristic(d, f, g);

        if (d1)
        {
           fmpz_poly_divrem_divconquer(q, r, d, a);

           result = fmpz_poly_is_zero(r) && _t_gcd_is_canonical(d);
           if (!result)
           {
              flint_printf("FAIL (check a | gcd(af, ag)):\n");
              flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n");
              flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n");
              flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n");
              flint_printf("d = "), fmpz_poly_print(d), flint_printf("\n");
              fflush(stdout);
              flint_abort();
           }
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    /* Check that a == GCD(af, ag) when GCD(f, g) = 1 */
    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, d, f, g, q, r;

        fmpz_poly_init(a);
        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 100) + 1, 200);
        do {
           fmpz_poly_randtest(f, state, n_randint(state, 100), 200);
           fmpz_poly_randtest(g, state, n_randint(state, 100), 200);
           fmpz_poly_gcd_heuristic(d, f, g);
        } while (!(d->length == 1 && fmpz_is_one(d->coeffs)));

        fmpz_poly_mul(f, a, f);
        fmpz_poly_mul(g, a, g);
        d1 = fmpz_poly_gcd_heuristic(d, f, g);

        if (d1)
        {
           if (!_t_gcd_is_canonical(a)) fmpz_poly_neg(a, a);

           result = fmpz_poly_equal(d, a) && _t_gcd_is_canonical(d);
           if (!result)
           {
              flint_printf("FAIL (check a == gcd(af, ag) when gcd(f, g) = 1):\n");
              flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n");
              flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n");
              flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n");
              flint_printf("d = "), fmpz_poly_print(d), flint_printf("\n");
              fflush(stdout);
              flint_abort();
           }
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    /*
	   Check that gcd(f, ga) divides f and ga for small generic f, g
	   and a small linear factor a. Exercises a bug found by Anton Mellit.
	*/
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, d, f, g, q, r;

		fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_init(a);
        fmpz_poly_randtest(f, state, n_randint(state, 10), 8);
        fmpz_poly_randtest(g, state, n_randint(state, 10), 4);

		/* multiply by small linear factor */
		fmpz_poly_set_coeff_si(a, 0, n_randint(state, 2) ? 1 : -1);
		fmpz_poly_set_coeff_si(a, 1, 1);
		fmpz_poly_mul(g, g, a);

        d1 = fmpz_poly_gcd_heuristic(d, f, g);

        if (d1)
        {
           if (fmpz_poly_is_zero(d))
		      result = fmpz_poly_is_zero(f) && fmpz_poly_is_zero(g);
		   else
		   {
		      fmpz_poly_divrem_divconquer(q, r, f, d);
              result = fmpz_poly_is_zero(r);
              fmpz_poly_divrem_divconquer(q, r, g, d);
              result &= fmpz_poly_is_zero(r);
		   }

           if (!result)
           {
              flint_printf("FAIL (gcd(f, g) | f and g):\n");
              flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n");
              flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n");
              flint_printf("d = "), fmpz_poly_print(d), flint_printf("\n");
              fflush(stdout);
              flint_abort();
           }
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    /* Sebastian's test case */
    {
       fmpz_poly_t a, b, d;

       fmpz_poly_init(a);
       fmpz_poly_init(b);
       fmpz_poly_init(d);

       fmpz_poly_set_coeff_ui(b, 2, 1);
       fmpz_poly_set_coeff_si(a, 0, -32);
       fmpz_poly_set_coeff_si(a, 1, 24);

       fmpz_poly_gcd_heuristic(d, a, b);

       result = (d->length == 1 && fmpz_is_one(d->coeffs));
       if (!result)
       {
          flint_printf("FAIL (check 1 == gcd(x^2, 24*x - 32):\n");
          fmpz_poly_print(d); flint_printf("\n");
          fflush(stdout);
          flint_abort();
       }

       fmpz_poly_clear(a);
       fmpz_poly_clear(b);
       fmpz_poly_clear(d);
    }

    /* Anton Mellit's test case */
    {
       fmpz_poly_t a, b, d;
       int heuristic;

       fmpz_poly_init(a);
       fmpz_poly_init(b);
       fmpz_poly_init(d);

	   /*
	       b = 3*q^12 - 8*q^11 - 24*q^10 - 48*q^9 - 84*q^8 - 92*q^7 - 92*q^6 -
               70*q^5 - 50*q^4 - 27*q^3 - 13*q^2 - 4*q - 1
		   a = q^13 - 2*q^12 + 2*q^10 - q^9
	   */
       fmpz_poly_set_str(b, "13  -1 -4 -13 -27 -50 -70 -92 -92 -84 -48 -24 -8 3");
	   fmpz_poly_set_str(a, "14  0 0 0 0 0 0 0 0 0 -1 2 0 -2 1");

       heuristic = fmpz_poly_gcd_heuristic(d, a, b);

       result = (heuristic == 0 || (d->length == 1 && fmpz_is_one(d->coeffs)));
       if (!result)
       {
          flint_printf("FAIL Mellit test case:\n");
          fmpz_poly_print(d); flint_printf("\n");
          fflush(stdout);
          flint_abort();
       }

       fmpz_poly_clear(a);
       fmpz_poly_clear(b);
       fmpz_poly_clear(d);
    }

    /* Daniel's test case */
  {
        fmpz_poly_t a, b, g;                                                                     fmpz_poly_init(a);                                                                       fmpz_poly_init(b);
        fmpz_poly_init(g);
        fmpz_poly_set_str(a, "40  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -7609399 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 44");
        fmpz_poly_set_str(b, "40  0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 54909036 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -59769402");
        fmpz_poly_gcd(g, a, b);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);
    }

  TEST_FUNCTION_END(state);
}
