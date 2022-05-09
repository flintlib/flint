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
   mp_limb_t d, r = 0;
   double dpre;
   ulong i;
   mp_ptr array = (mp_ptr) flint_malloc(1024*sizeof(mp_limb_t));
   FLINT_TEST_INIT(state);
   
   
   for (i = 0; i < count; i++)
   {
      int j;
      d = n_randtest(state);
      if (d == UWORD(0)) d++;

      dpre = n_precompute_inverse(d);

      for (j = 0; j < 1024; j++)
      {
         array[j] = n_randtest(state);
      }

      prof_start();
      for (j = 0; j < 10000; j++)
      {
         r += n_mod2_precomp(array[j&1023], d, dpre);  
      }
      prof_stop();
   }

   if (r == 0) flint_abort();

   flint_randclear(state);
   flint_free(array);
}

int main(void)
{
   double min, max;
   
   prof_repeat(&min, &max, sample, NULL);
   
   flint_printf("mod2_precomp min time is %.3f cycles, max time is %.3f cycles\n", 
           (min/(double)FLINT_CLOCK_SCALE_FACTOR)/10000.0, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/10000.0);

   return 0;
}
