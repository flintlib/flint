/*
    Copyright (C) 2015 Nitin Kumar
    Copyright (C) 2016 William Hart
    Copyright (C) 2020 Dan Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "qsieve.h"

void randprime(fmpz_t p, flint_rand_t state, slong bits)
{
    fmpz_randbits(p, state, bits);

    if (fmpz_sgn(p) < 0)
       fmpz_neg(p, p);

    if (fmpz_is_even(p))
       fmpz_add_ui(p, p, 1);

    while (!fmpz_is_probabprime(p))
       fmpz_add_ui(p, p, 2);
}

TEST_FUNCTION_START(qsieve_factor, state)
{
   slong i;
   fmpz_t n, x, y, z;
   fmpz_factor_t factors;
   slong max_threads = 5;
   slong tmul = 3;
#ifdef _WIN32
   tmul = 1;
#endif

   fmpz_init(x);
   fmpz_init(y);
   fmpz_init(z);
   fmpz_init(n);

   /* Test n with large prime factor */
   {
      fmpz_set_str(n, "12387192837918273918723981291837121933111751252512531193171", 10);

      fmpz_factor_init(factors);

      qsieve_factor(factors, n);

      if (factors->num < 5)
      {
         flint_printf("FAIL:\n");
         flint_printf("Test n with large prime factor\n");
         flint_printf("%ld factors found\n", factors->num);
         fflush(stdout);
         flint_abort();
      }

      fmpz_factor_clear(factors);
   }

   /* Test random n, two factors */
   for (i = 0; i < tmul*flint_test_multiplier(); i++)
   {
      slong bits = 40;

      randprime(x, state, bits);
      do {
         randprime(y, state, bits);
      } while (fmpz_equal(x, y));

      fmpz_mul(n, x, y);

      fmpz_factor_init(factors);

      flint_set_num_threads(n_randint(state, max_threads) + 1);

      qsieve_factor(factors, n);

      if (factors->num < 2)
      {
         flint_printf("FAIL:\n");
         flint_printf("Test random n, two factors\ni = %wd\n", i);
         flint_printf("%ld factors found\n", factors->num);
         fflush(stdout);
         flint_abort();
      }

      fmpz_factor_clear(factors);
   }

   /* Test random n, three factors */
   for (i = 0; i < tmul*flint_test_multiplier(); i++)
   {
      randprime(x, state, 40);
      do {
         randprime(y, state, 40);
      } while (fmpz_equal(x, y));
      do {
         randprime(z, state, 40);
      } while (fmpz_equal(x, z) || fmpz_equal(y, z));

      fmpz_mul(n, x, y);
      fmpz_mul(n, n, z);

      fmpz_factor_init(factors);

      flint_set_num_threads(n_randint(state, max_threads) + 1);

      qsieve_factor(factors, n);

      if (factors->num < 3)
      {
         flint_printf("FAIL:\n");
         flint_printf("Test random n, three factors\ni = %wd\n", i);
         flint_printf("%ld factors found\n", factors->num);
         fflush(stdout);
         flint_abort();
      }

      fmpz_factor_clear(factors);
   }

   /* Test random n, small factors */
   for (i = 0; i < tmul*flint_test_multiplier(); i++)
   {
      randprime(x, state, 10);
      do {
         randprime(y, state, 10);
      } while (fmpz_equal(x, y));
      randprime(z, state, 40);

      fmpz_mul(n, x, y);
      fmpz_mul(n, n, z);

      fmpz_factor_init(factors);

      flint_set_num_threads(n_randint(state, max_threads) + 1);

      qsieve_factor(factors, n);

      if (factors->num < 3)
      {
         flint_printf("FAIL:\n");
         flint_printf("Test random n, small factors\ni = %wd\n", i);
         flint_printf("%ld factors found\n", factors->num);
         fflush(stdout);
         flint_abort();
      }

      fmpz_factor_clear(factors);
   }

   fmpz_clear(n);
   fmpz_clear(x);
   fmpz_clear(y);
   fmpz_clear(z);

   TEST_FUNCTION_END(state);
}
