/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "test_helpers.h"

TEST_FUNCTION_START(smul_ppmm, state)
{
   int i, result;

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t ph1, pl1, ph2, pl2, pl2old, n1, n2, m1, m2, bit;
      int j, sign;

      n1 = n_randtest(state);
      n2 = n_randtest(state);

      smul_ppmm(ph1, pl1, n1, n2);

      m1 = n1;
      m2 = n2;

      sign = 1;
      if ((mp_limb_signed_t) m1 < WORD(0))
      {
         sign = -1;
         m1 = -m1;
      }

      if ((mp_limb_signed_t) m2 < WORD(0))
      {
         sign = -sign;
         m2 = -m2;
      }

      pl2old = UWORD(0);
      pl2 = UWORD(0);
      ph2 = UWORD(0);
      bit = UWORD(1);
      for (j = 0; j < FLINT_BITS; j++)
      {
         if (m2 & bit)
         {
            pl2 += (m1 << j);
            ph2 += (pl2 < pl2old);
            ph2 += r_shift(m1, FLINT_BITS - j);
            pl2old = pl2;
         }
         bit <<= 1;
      }

      if (sign == -1)
         sub_ddmmss(ph2, pl2, 0, 0, ph2, pl2);

      result = ((ph2 == ph1) && (pl2 == pl1));

      if (!result)
         flint_throw(FLINT_TEST_FAIL,
                 "Inputs: u = %wd, v = %wd\n\n"
                 "Got:      r1 = %wd, r0 = %wd\n"
                 "Expected: r1 = %wd, r0 = %wd\n",
                 n1, n2, ph1, pl1, ph2, pl2);
   }

   TEST_FUNCTION_END(state);
}
