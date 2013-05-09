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
} info_t;

void sample(void * arg, ulong count)
{
   mp_limb_t n;
   nmod_t mod;
   info_t * info = (info_t *) arg;
   mp_bitcnt_t bits = info->bits;
   mp_ptr vec = _nmod_vec_init(1000);
   mp_ptr vec2 = _nmod_vec_init(1000);
   mp_size_t j;
   len_t i;
   flint_rand_t state;
   flint_randinit(state);
    
   for (j = 0; j < 1000; j++)
      vec[j] = n_randlimb(state);

   prof_start();
   for (i = 0; i < count; i++)
   {
      n = n_randbits(state, bits);
      if (n == 0UL) n++;
      
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
   mp_bitcnt_t i;

   for (i = 2; i <= FLINT_BITS; i++)
   {
      info.bits = i;

	  prof_repeat(&min, &max, sample, (void *) &info);

      printf("bits %ld, c/l = %.1lf\n", 
         i, (min/(double)FLINT_CLOCK_SCALE_FACTOR)/1000
	  );
   }

   return 0;
}
