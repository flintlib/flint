/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
TEST_FUNCTION_START(fmpz_gcdinv, state)
{
   int i, result;

   /* Test a f == d mod g (specifically d > 1) */
   for (i = 0; i < 1000 * flint_test_multiplier(); i++)
   {
      fmpz_t d, a, f, g, t, x;
      int aliasing = n_randint(state, 5);

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

      if (aliasing == 0)
      {
          fmpz_gcdinv(d, a, f, g);
      }
      else if (aliasing == 1)
      {
          fmpz_set(d, f);
          fmpz_gcdinv(d, a, d, g);
      }
      else if (aliasing == 2)
      {
          fmpz_set(a, g);
          fmpz_gcdinv(d, a, f, a);
      }
      else if (aliasing == 3)
      {
          fmpz_set(a, f);
          fmpz_gcdinv(d, a, a, g);
      }
      else
      {
          fmpz_set(d, g);
          fmpz_gcdinv(d, a, f, d);
      }

      fmpz_mul(t, a, f);
      fmpz_mod(t, t, g);

      result = ((fmpz_equal(t, d) || fmpz_is_zero(f)) &&
	           fmpz_cmp_ui(a, 0) >= 0 && fmpz_cmp(a, g) < 0)
                && _fmpz_is_canonical(d) && _fmpz_is_canonical(a);
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

   TEST_FUNCTION_END(state);
}
