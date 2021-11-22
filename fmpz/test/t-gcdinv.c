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
#include "ulong_extras.h"
#include "fmpz.h"

int main(void)
{
   int i, result;
   FLINT_TEST_INIT(state);

   flint_printf("gcdinv....");
   fflush(stdout);

   /* Test aliasing of d and f, a and g */
   for (i = 0; i < 1000 * flint_test_multiplier(); i++)
   {
      fmpz_t d, a, f, g, F, G;

      fmpz_init(d);
      fmpz_init(a);
      fmpz_init(f);
      fmpz_init(g);
      fmpz_init(F);
      fmpz_init(G);

      fmpz_randtest_unsigned(G, state, 200);
      fmpz_add_ui(G, G, 1);
      fmpz_randm(F, state, G);
      fmpz_set(f, F);
      fmpz_set(g, G);

      fmpz_gcdinv(d, a, f, g);
      fmpz_gcdinv(f, g, f, g);

      result = (fmpz_equal(d, f) && (fmpz_equal(a, g) || fmpz_is_zero(F)));
      if (!result)
      {
         flint_printf("FAIL:\n\n");
         flint_printf("d = "), fmpz_print(d), flint_printf("\n");
         flint_printf("a = "), fmpz_print(a), flint_printf("\n");
         flint_printf("f = "), fmpz_print(f), flint_printf("\n");
         flint_printf("g = "), fmpz_print(g), flint_printf("\n");
         flint_printf("F = "), fmpz_print(F), flint_printf("\n");
         flint_printf("G = "), fmpz_print(G), flint_printf("\n");
         fflush(stdout);
         flint_abort();
      }

      fmpz_clear(d);
      fmpz_clear(a);
      fmpz_clear(f);
      fmpz_clear(g);
      fmpz_clear(F);
      fmpz_clear(G);
   }

   /* Test aliasing of a and f, d and g */
   for (i = 0; i < 1000 * flint_test_multiplier(); i++)
   {
      fmpz_t d, a, f, g, F, G;

      fmpz_init(d);
      fmpz_init(a);
      fmpz_init(f);
      fmpz_init(g);
      fmpz_init(F);
      fmpz_init(G);

      fmpz_randtest_unsigned(G, state, 200);
      fmpz_add_ui(G, G, 1);
      fmpz_randm(F, state, G);
      fmpz_set(f, F);
      fmpz_set(g, G);

      fmpz_gcdinv(d, a, f, g);
      fmpz_gcdinv(g, f, f, g);

      result = (fmpz_equal(d, g) && (fmpz_equal(a, f) || fmpz_is_zero(F)));
      if (!result)
      {
         flint_printf("FAIL:\n\n");
         flint_printf("d = "), fmpz_print(d), flint_printf("\n");
         flint_printf("a = "), fmpz_print(a), flint_printf("\n");
         flint_printf("f = "), fmpz_print(f), flint_printf("\n");
         flint_printf("g = "), fmpz_print(g), flint_printf("\n");
         flint_printf("F = "), fmpz_print(F), flint_printf("\n");
         flint_printf("G = "), fmpz_print(G), flint_printf("\n");
         fflush(stdout);
         flint_abort();
      }

      fmpz_clear(d);
      fmpz_clear(a);
      fmpz_clear(f);
      fmpz_clear(g);
      fmpz_clear(F);
      fmpz_clear(G);
   }

   /* Test a f == d mod g (generically d == 1) */
   for (i = 0; i < 1000 * flint_test_multiplier(); i++)
   {
      fmpz_t d, a, f, g, t;

      fmpz_init(d);
      fmpz_init(a);
      fmpz_init(f);
      fmpz_init(g);
      fmpz_init(t);

      fmpz_randtest_unsigned(g, state, 200);
      fmpz_add_ui(g, g, 1);
      fmpz_randm(f, state, g);

      fmpz_gcdinv(d, a, f, g);

      fmpz_mul(t, a, f);
      fmpz_mod(t, t, g);

      result = ((fmpz_equal(t, d) || fmpz_is_zero(f)) &&
		fmpz_cmp_ui(a, 0) >= 0 && fmpz_cmp(a, g) < 0);
      if (!result)
      {
         flint_printf("FAIL:\n\n");
         flint_printf("d = "), fmpz_print(d), flint_printf("\n");
         flint_printf("a = "), fmpz_print(a), flint_printf("\n");
         flint_printf("f = "), fmpz_print(f), flint_printf("\n");
         flint_printf("g = "), fmpz_print(g), flint_printf("\n");
         flint_printf("t = "), fmpz_print(t), flint_printf("\n");
         fflush(stdout);
         flint_abort();
      }

      fmpz_clear(d);
      fmpz_clear(a);
      fmpz_clear(f);
      fmpz_clear(g);
      fmpz_clear(t);
   }

   /* Test a f == d mod g (specifically d > 1) */
   for (i = 0; i < 1000 * flint_test_multiplier(); i++)
   {
      fmpz_t d, a, f, g, t, x;

      fmpz_init(d);
      fmpz_init(a);
      fmpz_init(f);
      fmpz_init(g);
      fmpz_init(t);
      fmpz_init(x);

      fmpz_randtest_unsigned(g, state, 200);
      fmpz_add_ui(g, g, 1);
      fmpz_randm(f, state, g);
      fmpz_randtest_unsigned(x, state, 100);
      fmpz_add_ui(x, x, 1);
      fmpz_mul(f, f, x);
      fmpz_mul(g, g, x);

      fmpz_gcdinv(d, a, f, g);

      fmpz_mul(t, a, f);
      fmpz_mod(t, t, g);

      result = ((fmpz_equal(t, d) || fmpz_is_zero(f)) &&
	           fmpz_cmp_ui(a, 0) >= 0 && fmpz_cmp(a, g) < 0);
      if (!result)
      {
         flint_printf("FAIL:\n\n");
         flint_printf("d = "), fmpz_print(d), flint_printf("\n");
         flint_printf("a = "), fmpz_print(a), flint_printf("\n");
         flint_printf("f = "), fmpz_print(f), flint_printf("\n");
         flint_printf("g = "), fmpz_print(g), flint_printf("\n");
         flint_printf("t = "), fmpz_print(t), flint_printf("\n");
         flint_printf("x = "), fmpz_print(x), flint_printf("\n");
         fflush(stdout);
         flint_abort();
      }

      fmpz_clear(d);
      fmpz_clear(a);
      fmpz_clear(f);
      fmpz_clear(g);
      fmpz_clear(t);
      fmpz_clear(x);
   }

   FLINT_TEST_CLEANUP(state);
    
   flint_printf("PASS\n");
   return 0;
}
