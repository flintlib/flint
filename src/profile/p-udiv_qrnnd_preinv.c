/*
    Copyright 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"

void sample(void * arg, ulong count)
{
   mp_limb_t d, q, r, dinv, norm;
   mp_ptr array = (mp_ptr) flint_malloc(200 * sizeof(mp_limb_t));
   ulong i;
   int j;

   FLINT_TEST_INIT(state);

   d = n_randtest_not_zero(state);
   norm = flint_clz(d);
   d <<= norm;

   for (i = 0; i < count; i++)
   {
      for (j = 0; j < 200; j+=2)
      {
         do
         {
            array[j] = n_randtest(state);
         } while (array[j] >= d);
         array[j + 1] = n_randtest(state);
      }

      dinv = n_preinvert_limb_prenorm(d);

      prof_start();
      for (j = 0; j < 200; j+=2)
      {
         udiv_qrnnd_preinv(q, r, array[j], array[j+1], d, dinv);
      }
      prof_stop();

      if (q + r == 0) flint_printf("\r");
   }

   flint_randclear(state);
   flint_free(array);
}

int main(void)
{
   double min, max;

   prof_repeat(&min, &max, sample, NULL);

   flint_printf("udiv_qrnnd_preinv min time is %.3f cycles, max time is %.3f cycles\n",
           (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

   return 0;
}
