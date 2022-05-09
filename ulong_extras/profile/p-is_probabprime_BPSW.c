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

typedef struct
{
   ulong bits;
} BPSW_t;

void sample(void * arg, ulong count)
{
   BPSW_t * params = (BPSW_t *) arg;
   ulong bits = params->bits;
   ulong i;
   mp_limb_t d;

   FLINT_TEST_INIT(state);
   

   for (i = 0; i < count; i++)
   {
      int j, res = 1;
      d = n_randbits(state, bits);
      while (!n_is_prime(d)) d++;

      prof_start();
      for (j = 0; j < 1000000; j++)
         res &= n_is_probabprime_BPSW(d);
      prof_stop();
      
      if (!res) flint_printf("Error\n");
   }

   flint_randclear(state);
}

int main(void)
{
   double min, max;
   BPSW_t params;
   int i;

   flint_printf("is_probabprime_BPSW:\n");
   
   for (i = 1; i <= 64; i++)
   {
      params.bits = i;
      prof_repeat(&min, &max, sample, &params);
      flint_printf("bits = %d, min time is %.3f cycles, max time is %.3f cycles\n", 
           i, (min/(double)FLINT_CLOCK_SCALE_FACTOR)/1000000, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/1000000);
   }

   return 0;
}
