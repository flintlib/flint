/*
    Copyright 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "ulong_extras.h"

void sample(void * arg, ulong count)
{
   ulong i;
   ulong a, b, d, r;
   double dpre;
   ulong_ptr array = (ulong_ptr) flint_malloc(1000*sizeof(ulong));
   FLINT_TEST_INIT(state);
   
   
   for (i = 0; i < count; i++)
   {
      int j;
      ulong bits = n_randint(state, 53) + 1;
      d = n_randbits(state, bits);
      a = n_randint(state, d);
      dpre = n_precompute_inverse(d);

      for (j = 0; j < 1000; j++)
      {
         array[i] = n_randint(state, d);
      }

      prof_start();
      for (j = 0; j < 1000; j++)
      {
         r = n_mulmod_precomp(a, array[j], d, dpre);
      }
      prof_stop();
   }

   flint_randclear(state);
   flint_free(array);
}

int main(void)
{
   double min, max;
   
   prof_repeat(&min, &max, sample, NULL);
   
   flint_printf("mulmod_precomp min time is %.3f cycles, max time is %.3f cycles\n", 
           (min/(double)FLINT_CLOCK_SCALE_FACTOR)/1000, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/1000);

   return 0;
}
