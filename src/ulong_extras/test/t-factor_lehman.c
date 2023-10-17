/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_factor_lehman, state)
{
   int i, result;

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random numbers */
   {
      mp_limb_t n1, n2;

      do
      {
         n1 = n_randtest_bits(state, n_randint(state, FLINT_BITS) + 1);
      } while (n_is_prime(n1) || (n1 < UWORD(2))
#if FLINT64 /* cannot compute enough primes */
         || (n1 >= UWORD(10000000000000000))
#endif
         );

      n2 = n_factor_lehman(n1);

      result = ((n1%n2) == UWORD(0) && n1 != n2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n1 = %wu, n2 = %wu\n", n1, n2);
         fflush(stdout);
         flint_abort();
      }
   }

   for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random products of two primes */
   {
      mp_limb_t n1, n2, n3, n, limit;

#if FLINT64
      limit = UWORD(100000000) - UWORD(100);
#else
      limit = UWORD(65535);
#endif

      n1 = n_randtest(state) % (limit + 1);
      n2 = n_randtest(state) % (limit + 1);

      n1 = n_nextprime(n1, 1);
      n2 = n_nextprime(n2, 1);

      /* test a specific bug */
#if FLINT64
      if (i == 0)
      {
           n1 = 72528697;
           n2 = 73339073;
      }
#endif

      n = n1*n2;

      n3 = n_factor_lehman(n);

      result = ((n%n3) == UWORD(0) && n != n3);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n = %wd, n3 = %wu\n", n, n3);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
