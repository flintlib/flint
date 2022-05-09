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

#define ITERS 1000

typedef struct
{
   ulong * composites;
   flint_bitcnt_t bits;
} fac_one_line_t;

void sample(void * arg, ulong count)
{
   fac_one_line_t * params = (fac_one_line_t *) arg;
   ulong i, j;
   mp_limb_t n2;
   
   for (i = 0; i < count; i++)
   {
      prof_start();
      for (j = 0; j < ITERS; j++)
	  {
		  /*n_factor_init(&factors);
	      n_factor(&factors, params->composites[j & 1023], 0);
		  if ((factors.num == 0) || (factors.num == 1 && factors.exp[0] < 2))
			  flint_printf("Error %wd\n", params->composites[j & 1023]);*/
        n2 = n_factor_lehman(params->composites[j & 1023]);
        if (n2 == params->composites[j & 1023])
           flint_printf("Error n = %wd\n", params->composites[j & 1023]);
	  }
	  prof_stop();
   }
}

void fill_array(ulong * ret, flint_bitcnt_t bits, flint_rand_t state)
{
   ulong n;
   n_factor_t factors;
   ulong i, primes = (1<<(bits/3))/10 + 1;
   
   for (i = 0; i < 1024; i++)
   {
	  do 
	  {
		 n_factor_init(&factors);
	     n = n_randbits(state, bits);
	  } while (n_is_probabprime(n) || (n_factor_trial(&factors, n, primes) != n));
	  ret[i] = n;
   }
      
}

int main(void)
{
   double min, max;
   fac_one_line_t params;
   FLINT_TEST_INIT(state);
   int i;
   

   params.composites = flint_malloc(1024*sizeof(ulong));

   flint_printf("factor_one_line:\n");
   
   for (i = 4; i <= 64; i++)
   {
      fill_array(params.composites, i, state);
      params.bits = i;
	  prof_repeat(&min, &max, sample, &params);
      flint_printf("bits = %d, time is %.3f us\n", 
		  i, max/(double)ITERS);
   }

   flint_randclear(state);
   flint_free(params.composites);
   return 0;
}
