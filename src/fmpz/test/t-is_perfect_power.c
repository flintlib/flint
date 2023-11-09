/*
    Copyright (C) 2009, 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_is_perfect_power, state)
{
   int i, result;
   ulong bits, exp;
   fmpz_t root, d, n, n2, pow;

   fmpz_init(d);
   fmpz_init(n);
   fmpz_init(n2);
   fmpz_init(pow);
   fmpz_init(root);

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that perfect powers pass the test */
   {
      bits = n_randint(state, 80) + 2;
      exp = n_randint(state, 10) + 2;
      fmpz_randtest(d, state, bits);
      fmpz_pow_ui(n, d, exp);

      result = fmpz_is_perfect_power(root, n);
      if (result == 0)
      {
         flint_printf("FAIL:\n");
         fmpz_print(n);
         flint_printf(" is declared not a perfect power\n");
         fflush(stdout);
         flint_abort();
      }

      fmpz_pow_ui(pow, root, result);
      if (!fmpz_equal(pow, n))
      {
         flint_printf("FAIL:\n");
         fmpz_print(root);
         flint_printf("^%d != ", result);
         fmpz_print(n);
         fflush(stdout);
         flint_abort();
      }
   }

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that non perfect powers fail */
   {
      mpz_t d_m;
      mpz_init(d_m);

      do
      {
         bits = n_randint(state, 1000) + 1;
         fmpz_randtest(n, state, bits);
         fmpz_get_mpz(d_m, n);
      } while (mpz_perfect_power_p(d_m));

      result = !fmpz_is_perfect_power(root, n);
      if (!result)
      {
         flint_printf("FAIL:\n");
         fmpz_print(n);
         flint_printf(" is declared a perfect power\n");
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(d_m);
   }

   for (i = 0; i < 100 * flint_test_multiplier(); i++) /* aliasing test, perfect powers */
   {
      bits = n_randint(state, 80) + 2;
      exp = n_randint(state, 10) + 2;
      fmpz_randtest(d, state, bits);
      fmpz_pow_ui(n, d, exp);
      fmpz_set(n2, n);

      result = fmpz_is_perfect_power(n, n);
      if (result == 0)
      {
         flint_printf("FAIL:\n");
         fmpz_print(n2);
         flint_printf(" is declared not a perfect power\n");
         fflush(stdout);
         flint_abort();
      }

      fmpz_pow_ui(pow, n, result);
      if (!fmpz_equal(pow, n2))
      {
         flint_printf("FAIL:\n");
         fmpz_print(n);
         flint_printf("^%d != ", result);
         fmpz_print(n2);
         fflush(stdout);
         flint_abort();
      }
   }

   for (i = 0; i < 100 * flint_test_multiplier(); i++) /* aliasing test, non perfect powers */
   {
      mpz_t d_m;
      mpz_init(d_m);

      do
      {
         bits = n_randint(state, 1000) + 1;
         fmpz_randtest(n, state, bits);
         fmpz_get_mpz(d_m, n);
      } while (mpz_perfect_power_p(d_m));

      result = !fmpz_is_perfect_power(n, n);
      if (!result)
      {
         flint_printf("FAIL:\n");
         fmpz_print(n2);
         flint_printf(" is declared a perfect power\n");
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(d_m);
   }

   fmpz_clear(n);
   fmpz_clear(n2);
   fmpz_clear(d);
   fmpz_clear(pow);
   fmpz_clear(root);

   TEST_FUNCTION_END(state);
}
