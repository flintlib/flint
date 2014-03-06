/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright 2014 Abhinav Baid

******************************************************************************/

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
   mp_bitcnt_t bits1;
   mp_bitcnt_t bits2;
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

void fill_array(ulong * ret, mp_bitcnt_t bits, flint_rand_t state)
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
   int i, j;
   

   params.rnums1 = flint_malloc(1024*sizeof(ulong));
   params.rnums2 = flint_malloc(1024*sizeof(ulong));

   flint_printf("n_gcd:\n");
   
   for (i = 1; i <= 64; i++)
   {
      fill_array(params.rnums1, i, state);
      params.bits1 = i;
      for(j = i + 1; j <= 64; j++) {
		  fill_array(params.rnums2, j, state);
		  prof_repeat(&min, &max, sample, &params);
		  flint_printf("bits1 = %d, bits2 = %d, time is %.3f us\n", 
						i, j, max/(double)ITERS);
	  }
	}

   flint_randclear(state);
   flint_free(params.rnums1);
   flint_free(params.rnums2);
   return 0;
}
