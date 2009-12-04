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

 (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"

typedef struct
{
   ulong bits;
} BPSW_t;

void sample(void * arg, ulong count)
{
   BPSW_t * params = (BPSW_t *) arg;

   ulong bits = params->bits;

   mp_limb_t n, d, r, norm;
   double dpre;
   
   for (ulong i = 0; i < count; i++)
   {
      d = n_randbits(bits);
      while (!n_is_prime(d)) d++;
      
      int res = 1;

      prof_start();
      for (mp_size_t i = 0; i < 1000000; i++)
         res &= n_is_probabprime_BPSW(d);
      prof_stop();
      
      if (!res) printf("Error\n");
   }
}

int main(void)
{
   double min, max;
   BPSW_t params;

   printf("is_probabprime_BPSW:\n");
   
   for (ulong i = 1; i <= 64; i++)
   {
      params.bits = i;
      prof_repeat(&min, &max, sample, &params);
      printf("bits = %ld, min time is %.3f cycles, max time is %.3f cycles\n", 
           i, (min/(double)FLINT_CLOCK_SCALE_FACTOR)/1000000, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/1000000);
   }

   return 0;
}
