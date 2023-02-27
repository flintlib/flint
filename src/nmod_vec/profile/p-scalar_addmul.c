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
   slong length;
} info_t;

void sample(void * arg, ulong count)
{
   mp_limb_t n, c;
   nmod_t mod;
   info_t * info = (info_t *) arg;
   flint_bitcnt_t bits = info->bits;
   slong length = info->length;
   slong i, j;
   mp_ptr vec = _nmod_vec_init(length);
   mp_ptr vec2 = _nmod_vec_init(length);
   FLINT_TEST_INIT(state);
   
    
   for (i = 0; i < count; i++)
   {
      n = n_randbits(state, bits);
      if (n == UWORD(0)) n++;
      c = n_randint(state, n);
      for (j = 0; j < length; j++)
         vec[j] = n_randint(state, n);
      for (j = 0; j < length; j++)
         vec2[j] = n_randint(state, n);

	  nmod_init(&mod, n);

      prof_start();
      for (j = 0; j < 30; j++)
		 _nmod_vec_scalar_addmul_nmod(vec2, vec, length, c, mod);
	  prof_stop();
   }
   
   flint_randclear(state);
   _nmod_vec_clear(vec);
   _nmod_vec_clear(vec2);
}

int main(void)
{
   double min0, min1, min2, max;
   info_t info;
   flint_bitcnt_t i;

   for (i = 2; i <= FLINT_BITS; i++)
   {
      info.bits = i;

	  info.length = 4;
	  prof_repeat(&min0, &max, sample, (void *) &info);

	  info.length = 1024;
	  prof_repeat(&min1, &max, sample, (void *) &info);

	  info.length = 65536;
	  prof_repeat(&min2, &max, sample, (void *) &info);

      flint_printf("bits %wd, length 4 %.1lf c/l, length 128 %.1lf c/l, length 65536 %.1lf c/l\n", 
         i,
         (min0/(double)FLINT_CLOCK_SCALE_FACTOR)/(4*30),
         (min1/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*30),
		 (min2/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*30)
	  );
   }

   return 0;
}
