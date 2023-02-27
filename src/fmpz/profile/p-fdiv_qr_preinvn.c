/*
    Copyright 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

typedef struct
{
   slong limbs;
   int algo;
} info_t;

void sample(void * arg, ulong count)
{
   info_t * info = (info_t *) arg;
   slong limbs = info->limbs, i, j;
   int algo = info->algo;
   int scale = 200;

   FLINT_TEST_INIT(state);
   

   fmpz_t a, b, c, r;
   fmpz_preinvn_t inv;
   
   fmpz_init(a);
   fmpz_init(b);
   fmpz_init(c);
   fmpz_init(r);
           
   for (i = 0; i < count; i++)
   {
      fmpz_randbits(a, state, (2*limbs - 1)*FLINT_BITS);
      fmpz_randbits(b, state, limbs*FLINT_BITS);
      
      fmpz_preinvn_init(inv, b);
	
      prof_start();
      if (algo == 1)
      {
         for (j = 0; j < scale; j++)
         {
            fmpz_fdiv_qr_preinvn(c, r, a, b, inv);
         }
      } else
      {
         for (j = 0; j < scale; j++)
         {
            fmpz_fdiv_qr(c, r, a, b);
         }
     }
	   prof_stop();
   }
  
   fmpz_preinvn_clear(inv);
   fmpz_clear(a);
   fmpz_clear(b);
   fmpz_clear(c);
   fmpz_clear(r);
   flint_randclear(state);
}

int main(void)
{
   double min, max;
   info_t info;
   slong k, scale;

   printf("1: With precomputed inverse\n");
   printf("2: Without precomputed inverse\n\n");

   for (k = 1; k <= 10000; k = (slong) ceil(1.1*k))
   {
      info.limbs = k;
      info.algo = 1;

      scale = 200;
   
      prof_repeat(&min, &max, sample, (void *) &info);
         
      flint_printf("1: limbs %wd, min %.3g ms, max %.3g ms\n", 
           info.limbs,
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0
	     );

     info.algo = 2;
     
     prof_repeat(&min, &max, sample, (void *) &info);
         
      flint_printf("2: limbs %wd, min %.3g ms, max %.3g ms\n\n", 
           info.limbs,
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0
	     );
   }

   return 0;
}
