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
   ulong i;
   mp_ptr arr = (mp_ptr) flint_malloc(1024*sizeof(mp_limb_t));
   FLINT_TEST_INIT(state);
   
      
   for (i = 0; i < count; i++)
   {
      int j;
      d = n_randbits(state, bits);
      if (d == UWORD(0)) d++;

      dinv = n_preinvert_limb(d);
      dpre = n_precompute_inverse(d);

      for (j = 0; j < 1024; j++)
      {
         arr[j] = n_randbits(state, FLINT_BITS);
      }

	  switch (type)
	  {
	  /*case 1:

         prof_start();
         for (mp_size_t j = 0; j < UWORD(10000); j++)
         {
            r += n_empty(arr[j&1023], d, dinv);  
         }
	     prof_stop();

	  break;*/

	  case 2:

         prof_start();
         for (j = 0; j < 10000; j++)
         {
            r += n_mod2_preinv(arr[j&1023], d, dinv);  
         }
	     prof_stop();

	  break;

	  /*case 3:

         prof_start();
         for (mp_size_t j = 0; j < UWORD(10000); j++)
         {
            r += n_mod3_preinv(arr[j&1023], d, dinv);  
         }
	     prof_stop();

	  break;*/

	  case 4:

         prof_start();
         for (j = 0; j < 10000; j++)
         {
            r += n_mod2_precomp(arr[j&1023], d, dpre);  
         }
	     prof_stop();

	  break;

	  case 5:

         prof_start();
         for (j = 0; j < 10000; j++)
         {
            r += n_mod_precomp(arr[j&1023], d, dpre);  
         }
	     prof_stop();

	  break;
	  }

   }
  
   if (r == UWORD(9879875897)) abort();

   flint_randclear(state);
   flint_free(arr);
}

int main(void)
{
   double min1, min2, min3, min4, min5, max;
   info_t info;
   int i;

   for (i = 1; i <= FLINT_BITS; i++)
   {
      info.bits = i;

      /*
	     info.type = 1;
         prof_repeat(&min1, &max, sample, (void *) &info);
       */

	  info.type = 2;
      prof_repeat(&min2, &max, sample, (void *) &info);

      /*
	     info.type = 3;
         prof_repeat(&min3, &max, sample, (void *) &info);
       */

	  info.type = 4;
      prof_repeat(&min4, &max, sample, (void *) &info);

	  if (i >= 32 && i <= 53)
	  {
	     info.type = 5;
         prof_repeat(&min5, &max, sample, (void *) &info);

         flint_printf("bits %d, inv2 %.1f c/l, pre2 %.1f c/l, pre %.1f c/l\n", 
           i,
		   /* (min1/(double)FLINT_CLOCK_SCALE_FACTOR)/10000, */
           (min2/(double)FLINT_CLOCK_SCALE_FACTOR)/10000,
           /* (min3/(double)FLINT_CLOCK_SCALE_FACTOR)/10000, */
           (min4/(double)FLINT_CLOCK_SCALE_FACTOR)/10000,
           (min5/(double)FLINT_CLOCK_SCALE_FACTOR)/10000
	     );

	  } else
	  {
         flint_printf("bits %d, inv2 %.1f c/l, pre2 %.1f c/l\n", 
           i,
		   /* (min1/(double)FLINT_CLOCK_SCALE_FACTOR)/10000, */
           (min2/(double)FLINT_CLOCK_SCALE_FACTOR)/10000,
           /* (min3/(double)FLINT_CLOCK_SCALE_FACTOR)/10000, */
           (min4/(double)FLINT_CLOCK_SCALE_FACTOR)/10000
	     );
	  }
   }

   return 0;
}
