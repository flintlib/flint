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
#include <mpir.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"

void sample(void * arg, ulong count)
{
   mp_limb_t d;
   
   mp_limb_t * array = (mp_limb_t *) malloc(200*sizeof(mp_limb_t));

   d = n_randtest_not_zero();
      
   for (ulong i = 0; i < count; i++)
   {
      for (ulong i = 0; i < 200; i+=2)
      {
         do
         {
            array[i] = n_randtest();
         } while (array[i] >= d);
         array[i + 1] = n_randtest();  
      }
       
      prof_start();
      for (ulong i = 0; i < 200; i+=2)
      {
         udiv_qrnnd(array[i], array[i+1], array[i], array[i+1], d);
      }
      prof_stop();

      for (ulong i = 0; i < 200; i++)
         if (array[i] == 0) printf("\r");
   }

   free(array);
}

int main(void)
{
   double min, max;
   
   prof_repeat(&min, &max, sample, NULL);
   
   printf("udiv_qrnnd min time is %.3f cycles, max time is %.3f cycles\n", 
           (min/(double)FLINT_CLOCK_SCALE_FACTOR)/100, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/100);

   return 0;
}
