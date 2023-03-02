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
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

typedef struct
{
   flint_bitcnt_t bits;
} info_t;

void sample(void * arg, ulong count)
{
   mp_limb_t n;
   nmod_t mod;
   info_t * info = (info_t *) arg;
   flint_bitcnt_t bits = info->bits;
   mp_ptr vec = _nmod_vec_init(1000);
   mp_ptr vec2 = _nmod_vec_init(1000);
   mp_size_t j;
   slong i;
   FLINT_TEST_INIT(state);
   
    
   for (j = 0; j < 1000; j++)
      vec[j] = n_randlimb(state);

   prof_start();
   for (i = 0; i < count; i++)
   {
      n = n_randbits(state, bits);
      if (n == UWORD(0)) n++;
      
	  nmod_init(&mod, n);
      _nmod_vec_reduce(vec2, vec, 1000, mod);
   }
   prof_stop();
 
   flint_randclear(state);
   _nmod_vec_clear(vec);
   _nmod_vec_clear(vec2);
}

int main(void)
{
   double min, max;
   info_t info;
   flint_bitcnt_t i;

   for (i = 2; i <= FLINT_BITS; i++)
   {
      info.bits = i;

	  prof_repeat(&min, &max, sample, (void *) &info);

      flint_printf("bits %wd, c/l = %.1lf\n", 
         i, (min/(double)FLINT_CLOCK_SCALE_FACTOR)/1000
	  );
   }

   return 0;
}
