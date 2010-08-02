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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int i, result;
   ulong count = 0UL;
   printf("factor_SQUFOF....");
   fflush(stdout);

   for (i = 0; i < 3000; i++) /* Test random numbers */
   {
      mp_limb_t n1, n2, bits;

      do
         n1 = n_randtest();
      while (n_is_prime(n1) || (n1 == 1UL));

#if FLINT64
      n2 = n_factor_SQUFOF(n1, 10000);
#else
      n2 = n_factor_SQUFOF(n1, 2000);
#endif
      
      if (n2)
      {
         count++;
         result = ((n1%n2) == 0UL);

         if (!result)
         {
            printf("FAIL:\n");
            printf("n1 = %lu, n2 = %lu\n", n1, n2); 
            abort();
         }
      }
   }
   
   if (count < 2800)
   {
      printf("FAIL:\n");
      printf("Only %lu of 10000 numbers factored\n", count); 
      abort();
   }

   printf("PASS\n");
   return 0;
}
