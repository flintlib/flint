/*============================================================================

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

===============================================================================*/
/******************************************************************************

 (C) 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"

typedef struct
{
   mp_bitcnt_t bits;
   ulong type;
} info_t;

void sample(void * arg, ulong count)
{
   mp_limb_t n, d, dinv, r = 0, norm;
   double dpre;

   info_t * info = (info_t *) arg;

   mp_bitcnt_t bits = info->bits;
   ulong type = info->type;

   mp_limb_t * arr = malloc(1024*sizeof(mp_limb_t));
   mp_limb_t * arr2 = malloc(1024*sizeof(mp_limb_t));
      
   for (ulong i = 0; i < count; i++)
   {
      d = n_randbits(bits);
      if (d == 0UL) d++;

      dinv = n_preinvert_limb(d);
      
      for (mp_size_t j = 0; j < 1024; j++)
      {
         arr[j] = n_randbits(FLINT_BITS);
         arr2[j] = n_randint(n);
      }

	  switch (type)
	  {
	  case 1:

         prof_start();
         for (mp_size_t j = 0; j < 10000UL; j++)
         {
            r += n_lll_mod_preinv(arr2[j&1023], arr[j&1023], arr[(j+1)&1023], d, dinv);  
         }
	     prof_stop();

	  break;
	  }

   }
  
   if (r == 9879875897UL) abort();

   free(arr);
   free(arr2);
}

int main(void)
{
   double min1, min2, min3, min4, min5, max;
   
   info_t info;

   for (ulong i = FLINT_BITS/2 + 1; i <= FLINT_BITS; i++)
   {
      info.bits = i;

	  info.type = 1;
      prof_repeat(&min1, &max, sample, (void *) &info);

	  printf("bits %ld, ll_inv %.1f c/l\n", 
           i,
		   (min1/(double)FLINT_CLOCK_SCALE_FACTOR)/10000
 	  );
   }

   return 0;
}
