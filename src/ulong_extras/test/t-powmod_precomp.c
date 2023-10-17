/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_powmod_precomp, state)
{
   int i, result;

   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t a, d, r1, r2, bits;
      mpz_t a_m, d_m, r2_m;
      mp_limb_signed_t exp;
      double dpre;

      mpz_init(a_m);
      mpz_init(d_m);
      mpz_init(r2_m);

      bits = n_randint(state, FLINT_D_BITS) + 1;
      d = n_randtest_bits(state, bits);
      do
      {
         a = n_randtest(state) % d;
      } while (n_gcd(d, a) != UWORD(1));
      exp = n_randtest(state);

      dpre = n_precompute_inverse(d);
      r1 = n_powmod_precomp(a, exp, d, dpre);

      flint_mpz_set_ui(a_m, a);
      flint_mpz_set_ui(d_m, d);
      if (exp < WORD(0))
      {
         flint_mpz_powm_ui(r2_m, a_m, -exp, d_m);
         mpz_invert(r2_m, r2_m, d_m);
      } else
         flint_mpz_powm_ui(r2_m, a_m, exp, d_m);
      r2 = flint_mpz_get_ui(r2_m);

      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, exp = %wd, d = %wu\n", a, exp, d);
         flint_printf("r1 = %wu, r2 = %wu\n", r1, r2);
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(a_m);
      mpz_clear(d_m);
      mpz_clear(r2_m);
   }

   /* check 0^0 = 1 */
   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t bits, d, r;
      double dpre;

      bits = n_randint(state, FLINT_D_BITS) + 1;
      d = n_randtest_bits(state, bits);
      if (d == 0) d++;

      dpre = n_precompute_inverse(d);
      r = n_powmod_precomp(0, 0, d, dpre);

      result = (r == 1 || (d == 1 && r == 0));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("0^0 != 1 mod %wd\n", d);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
