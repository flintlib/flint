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

    Copyright 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"

void sample(void * arg, ulong count)
{
   mp_limb_t a, b, d, r, norm, dinv;
   mp_ptr array = (mp_ptr) flint_malloc(1000*sizeof(mp_limb_t));
   ulong i;
   FLINT_TEST_INIT(state);
   
   
   for (i = 0; i < count; i++)
   {
      int j;
      mp_limb_t bits = n_randint(state, 53) + 1;
      d = n_randbits(state, bits);
      a = n_randint(state, d);
      dinv = n_preinvert_limb(d);

      for (j = 0; j < 1000; j++)
      {
         array[j] = n_randint(state, d);
      }

      prof_start();
      for (j = 0; j < 1000; j++)
      {
         r = n_mulmod2_preinv(a, array[j], d, dinv);
      }
      prof_stop();
   }

   flint_randclear(state);
   flint_free(array);
}

int main(void)
{
   double min, max;
   
   prof_repeat(&min, &max, sample, NULL);
   
   flint_printf("mulmod2_precomp min time is %.3f cycles, max time is %.3f cycles\n", 
           (min/(double)FLINT_CLOCK_SCALE_FACTOR)/1000, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/1000);

   return 0;
}
