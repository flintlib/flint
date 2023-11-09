/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_factor.h"

TEST_FUNCTION_START(fmpz_factor_pp1, state)
{
   int i, j, result;
   ulong count = UWORD(0);
   gmp_randstate_t st;

   gmp_randinit_default(st);

   for (i = 0; i < 50 * flint_test_multiplier(); i++) /* Test random numbers */
   {
      flint_bitcnt_t bits;
      mpz_t m, n;
      fmpz_t n1, n2, r;

      mpz_init(n);
      mpz_init(m);
      fmpz_init(n1);
      fmpz_init(n2);
      fmpz_init(r);

      do {
         mpz_urandomb(n, st, n_randint(state, 128) + 2);
      } while (flint_mpz_cmp_ui(n, 2) < 0);
      do {
         mpz_urandomb(m, st, n_randint(state, 50) + 2);
      } while (!mpz_probab_prime_p(m, 20));
      mpz_mul(n, n, m);

      fmpz_set_mpz(n1, n);
      bits = FLINT_MIN(fmpz_bits(n1), FLINT_BITS);

      for (j = 0; j < 20; j++)
      {
         fmpz_factor_pp1(n2, n1, 10000, 10000, n_randbits(state, bits - 2) + 3);
         if (fmpz_cmp_ui(n2, 1) > 0) break;
      }

      if (fmpz_cmp_ui(n2, 1) > 0)
      {
         count++;
         fmpz_mod(r, n1, n2);
         result = (fmpz_is_zero(r));
         if (!result)
         {
            flint_printf("FAIL:\n");
            flint_printf("n1 = ");
            fmpz_print(n1);
            flint_printf(", n2 = ");
            fmpz_print(n2);
            flint_printf("\n");
            fmpz_print(r); flint_printf("\n");
            fflush(stdout);
            flint_abort();
         }
      }

      fmpz_clear(n1);
      fmpz_clear(n2);
      fmpz_clear(r);
      mpz_clear(m);
      mpz_clear(n);
   }

   if (count < 49 * flint_test_multiplier())
   {
      flint_printf("FAIL:\n");
      flint_printf("Only %wu numbers factored\n", count);
      fflush(stdout);
      flint_abort();
   }

   gmp_randclear(st);

   TEST_FUNCTION_END(state);
}
