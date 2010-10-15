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
/****************************************************************************

   Copyright (C) 2010 William Hart
   Copyright (C) 2010 Daniel Woodhouse

*****************************************************************************/

#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_mpoly.h"
#include "ulong_extras.h"

int main(void)
{
   long i;
   
   printf("init/init2/realloc/clear....");
   fflush(stdout);
   
   for (i = 0; i < 10000UL; i++) 
   {
      nmod_mpoly_t a;

	   mp_limb_t n = n_randtest_not_zero();

      nmod_mpoly_init2(a, n, n_randint(100), (long) 5, (ulong) 5);
      nmod_mpoly_clear(a);      
   }

   for (i = 0; i < 10000UL; i++) 
   {
      nmod_mpoly_t a;
      

	   mp_limb_t n = n_randtest_not_zero();

      nmod_mpoly_init2(a, n, n_randint(100), (long) 5, (ulong) 5);
      nmod_mpoly_realloc(a, n_randint(100));
      nmod_mpoly_realloc(a, n_randint(100));
      nmod_mpoly_clear(a);      
   }
   
   for (i = 0; i < 10000UL; i++) 
   {
      nmod_mpoly_t a;

	   mp_limb_t n = n_randtest_not_zero();

      nmod_mpoly_init(a, n, (long) 5, (ulong) 5);
      nmod_mpoly_randtest(a, n_randint(100));
      
      nmod_mpoly_clear(a);
   }
   
   printf("PASS\n");
   return 0;
}
