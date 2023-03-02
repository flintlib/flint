/*
    Copyright 2014 Abhinav Baid

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

#define ITERS 1000

typedef struct
{
   ulong * rnums1;
   ulong * rnums2;
   flint_bitcnt_t bits1;
   flint_bitcnt_t bits2;
} gcd_t;

void sample(void * arg, ulong count)
{
   gcd_t * params = (gcd_t *) arg;
   ulong i, j;
   
   for (i = 0; i < count; i++)
   {
      prof_start();
      for (j = 0; j < ITERS; j++)
	  {
        n_gcd(params->rnums2[j & 1023], params->rnums1[j & 1023]);
	  }
	  prof_stop();
   }
}

void fill_array(ulong * ret, flint_bitcnt_t bits, flint_rand_t state)
{
   ulong n;
   ulong i;
   
   for (i = 0; i < 1024; i++)
   {
	  n = n_randbits(state, bits);
	  ret[i] = n;
   }
}

int main(void)
{
   double min, max;
   gcd_t params;
   FLINT_TEST_INIT(state);
   int i;
   

   params.rnums1 = flint_malloc(1024*sizeof(ulong));
   params.rnums2 = flint_malloc(1024*sizeof(ulong));

   flint_printf("n_gcd:\n");
   
   for (i = 1; i <= 64; i++)
   {
      fill_array(params.rnums1, i, state);
      params.bits1 = i;
		  fill_array(params.rnums2, i, state);
		  prof_repeat(&min, &max, sample, &params);
		  flint_printf("bits1 = %d, bits2 = %d, time is %.3f us\n", 
						i, i, max/(double)ITERS);
	}

   flint_randclear(state);
   flint_free(params.rnums1);
   flint_free(params.rnums2);
   return 0;
}
