/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_is_perfect_power235, state)
{
   int i, result;
   ulong bits;
   ulong d;

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that square pass the test */
   {
       bits = n_randint(state, FLINT_BITS/2) + 1;
       d = n_randtest_bits(state, bits);

       result = n_is_perfect_power235(n_pow(d, 2));
       if (!result)
           TEST_FUNCTION_FAIL("d^2 = %wu is declared not a perfect power\n", d * d);
   }

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that cubes pass the test */
   {
       bits = n_randint(state, FLINT_BITS/3) + 1;
       d = n_randtest_bits(state, bits);

       result = n_is_perfect_power235(n_pow(d, 3));

       if (!result)
           TEST_FUNCTION_FAIL("d^3 = %wu is declared not a perfect power\n", d * d * d);
   }

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that fifth powers pass the test */
   {
       bits = n_randint(state, FLINT_BITS/5) + 1;
       d = n_randtest_bits(state, bits);

       result = n_is_perfect_power235(n_pow(d, 5));

       if (!result)
           TEST_FUNCTION_FAIL("d^5 = %wu is declared not a perfect power\n", d * d * d * d * d);
   }

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that non prefect powers fail */
   {
       mpz_t d_m;
       mpz_init(d_m);

       do
       {
           d = n_randtest(state);
           flint_mpz_set_ui(d_m, d);
       } while (mpz_perfect_power_p(d_m));

       result = !n_is_perfect_power235(d);
       if (!result)
           flint_printf("d = %wu is declared a perfect power\n", d);

       mpz_clear(d_m);
   }

   TEST_FUNCTION_END(state);
}
