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

typedef struct
{
   flint_bitcnt_t bits;
   ulong type;
} info_t;

void sample(void * arg, ulong count)
{
   mp_limb_t d, dinv, r = 0;
   info_t * info = (info_t *) arg;
   flint_bitcnt_t bits = info->bits;
   ulong type = info->type;
   ulong i;
   FLINT_TEST_INIT(state);
   
      
   mp_ptr arr  = (mp_ptr) flint_malloc(1024*sizeof(mp_limb_t));
   mp_ptr arr2 = (mp_ptr) flint_malloc(1024*sizeof(mp_limb_t));
      
   for (i = 0; i < count; i++)
   {
      int j;
      d = n_randbits(state, bits);
      if (d == UWORD(0)) d++;

      dinv = n_preinvert_limb(d);
      
      for (j = 0; j < 1024; j++)
      {
         arr[j] = n_randbits(state, FLINT_BITS);
         arr2[j] = n_randint(state, d);
      }

	  switch (type)
	  {
	  case 1:

         prof_start();
         for (mp_size_t j = 0; j < UWORD(10000); j++)
         {
            r += n_lll_mod_preinv(arr2[j&1023], arr[j&1023], arr[(j+1)&1023], d, dinv);  
         }
	     prof_stop();

	  break;
	  }

   }
  
   if (r == UWORD(9879875897)) flint_abort();

   flint_randclear(state);
   flint_free(arr);
   flint_free(arr2);
}

int main(void)
{
   double min1, max;
   info_t info;
   int i;

   for (i = FLINT_BITS/2 + 1; i <= FLINT_BITS; i++)
   {
      info.bits = i;
	  info.type = 1;
      prof_repeat(&min1, &max, sample, (void *) &info);

	  flint_printf("bits %d, ll_inv %.1f c/l\n", 
           i,
		   (min1/(double)FLINT_CLOCK_SCALE_FACTOR)/10000
 	  );
   }

   return 0;
}
