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

    Copyright 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

typedef struct
{
   mp_bitcnt_t bits;
   len_t length;
} info_t;

void sample(void * arg, ulong count)
{
   mp_limb_t n, c;
   nmod_t mod;
   info_t * info = (info_t *) arg;
   mp_bitcnt_t bits = info->bits;
   len_t length = info->length;
   len_t i, j;
   mp_ptr vec = _nmod_vec_init(length);
   mp_ptr vec2 = _nmod_vec_init(length);
   flint_rand_t state;
   flint_randinit(state);
    
   for (i = 0; i < count; i++)
   {
      n = n_randbits(state, bits);
      if (n == 0UL) n++;
      c = n_randint(state, n);
      for (j = 0; j < length; j++)
         vec[j] = n_randint(state, n);
      
	  nmod_init(&mod, n);

      prof_start();
      for (j = 0; j < 30; j++)
		 _nmod_vec_scalar_mul_nmod(vec2, vec, length, c, mod);
	  prof_stop();
   }
   
   flint_randclear(state);
   _nmod_vec_clear(vec);
   _nmod_vec_clear(vec2);
}

int main(void)
{
   double min1, min2, max;
   info_t info;
   mp_bitcnt_t i;

   for (i = 2; i <= FLINT_BITS; i++)
   {
      info.bits = i;

	  info.length = 1024;
	  prof_repeat(&min1, &max, sample, (void *) &info);

	  info.length = 65536;
	  prof_repeat(&min2, &max, sample, (void *) &info);

      printf("bits %ld, length 128 %.1lf c/l, length 65536 %.1lf c/l\n", 
         i, (min1/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*30),
		 (min2/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*30)
	  );
   }

   return 0;
}
