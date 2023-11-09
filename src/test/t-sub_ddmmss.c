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

TEST_FUNCTION_START(sub_ddmmss, state)
{
   int i, result;

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t dh1, dl1, dh2, dl2, mh, ml, sh, sl;

      mh = n_randtest(state);
      ml = n_randtest(state);
      sh = n_randtest(state);
      sl = n_randtest(state);

      if (n_randint(state, 10) == 0)
          sub_ddmmss(dh1, dl1, (slong) mh, (slong) ml, (slong) sh, (slong) sl);
      else
          sub_ddmmss(dh1, dl1, mh, ml, sh, sl);

      dl2 = ml - sl;
      dh2 = -(sl > ml);
      dh2 += mh;
      dh2 -= sh;

      result = ((dh2 == dh1) && (dl2 == dl1));

      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("mh = %wu, ml = %wu, sh = %wu, sl = %wu\n", mh, ml, sh, sl);
         flint_printf("dh2 = %wu, dh1 = %wu, dl2 = %wu, dl1 = %wu\n", dh2, dh1, dl2, dl1);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
