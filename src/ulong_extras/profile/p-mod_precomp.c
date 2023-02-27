/*
    Copyright 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"

void sample(void * arg, ulong count)
{
   mp_limb_t d, bits;
   double dpre;
   ulong i;
   mp_ptr array = (mp_ptr) flint_malloc(1000*sizeof(mp_limb_t));
   FLINT_TEST_INIT(state);
   
   
   for (i = 0; i < count; i++)
   {
      int j;
      bits = n_randint(state, 53) + 1;
      d = n_randbits(state, bits);
      dpre = n_precompute_inverse(d);

      for (j = 0; j < 1000; j++)
      {
         if (bits <= 32) array[j] = n_randint(state, d*d);
         else array[j] = n_randtest(state);
      }

      prof_start();
      for (j = 0; j < 1000; j++)
      {
         array[j] = n_mod_precomp(array[j], d, dpre);  
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
   
   flint_printf("mod_precomp min time is %.3f cycles, max time is %.3f cycles\n", 
           (min/(double)FLINT_CLOCK_SCALE_FACTOR)/1000, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/1000);

   return 0;
}
