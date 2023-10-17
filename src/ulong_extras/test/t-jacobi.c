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

TEST_FUNCTION_START(n_jacobi, state)
{
   int i, result;

   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t d;
      mpz_t a_m, d_m;
      mp_limb_signed_t a;
      int r1, r2;

      mpz_init(a_m);
      mpz_init(d_m);

      a = n_randtest(state);
      d = n_randtest_not_zero(state) | WORD(1);

      r1 = n_jacobi(a, d);

      flint_mpz_set_si(a_m, a);
      flint_mpz_set_ui(d_m, d);
      r2 = mpz_jacobi(a_m, d_m);

      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, d = %wu\n", a, d);
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(a_m);
      mpz_clear(d_m);
   }

   TEST_FUNCTION_END(state);
}
