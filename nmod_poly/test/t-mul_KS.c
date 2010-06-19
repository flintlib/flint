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

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int main(void)
{
   int result;
   printf("mul_KS....");
   fflush(stdout);
   
   // check aliasing of a and b
   for (ulong i = 0; i < 2000UL; i++) 
   {
      nmod_poly_t a, b, c;
      
	  mp_limb_t n = n_randtest_not_zero();
      
      nmod_poly_init(a, n);
      nmod_poly_init(b, n);
      nmod_poly_init(c, n);
      nmod_poly_randtest(b, n_randint(50));
      nmod_poly_randtest(c, n_randint(50));
   
	  nmod_poly_mul_KS(a, b, c, 0);
      nmod_poly_mul_KS(b, b, c, 0);
      
      result = (nmod_poly_equal(a, b));
      if (!result)
      {
         printf("Error:\n");
         nmod_poly_print(a); printf("\n\n");
         nmod_poly_print(b); printf("\n\n");
         abort();
      }

      nmod_poly_clear(a);
      nmod_poly_clear(b);
      nmod_poly_clear(c);
   }
   
   // check aliasing of a and c
   for (ulong i = 0; i < 2000UL; i++) 
   {
      nmod_poly_t a, b, c;
      
	  mp_limb_t n = n_randtest_not_zero();
      
      nmod_poly_init(a, n);
      nmod_poly_init(b, n);
      nmod_poly_init(c, n);
      nmod_poly_randtest(b, n_randint(50));
      nmod_poly_randtest(c, n_randint(50));
   
	  nmod_poly_mul_KS(a, b, c, 0);
      nmod_poly_mul_KS(c, b, c, 0);
      
      result = (nmod_poly_equal(a, c));
      if (!result)
      {
         printf("Error:\n");
         nmod_poly_print(a); printf("\n\n");
         nmod_poly_print(c); printf("\n\n");
         abort();
      }

      nmod_poly_clear(a);
      nmod_poly_clear(b);
      nmod_poly_clear(c);
   }
   
   // compare with mul_classical
   for (ulong i = 0; i < 2000UL; i++) 
   {
      nmod_poly_t a1, a2, b, c;
      
	  mp_limb_t n = n_randtest_not_zero();
      
      nmod_poly_init(a1, n);
      nmod_poly_init(a2, n);
      nmod_poly_init(b, n);
      nmod_poly_init(c, n);
      nmod_poly_randtest(b, n_randint(50));
      nmod_poly_randtest(c, n_randint(50));
      
	  nmod_poly_mul_classical(a1, b, c);
	  nmod_poly_mul_KS(a2, b, c, 0);
       
      result = (nmod_poly_equal(a1, a2));
      if (!result)
      {
		 printf("Error:\n");
		 nmod_poly_print(a1); printf("\n\n");
         nmod_poly_print(a2); printf("\n\n");
         abort();
      }

      nmod_poly_clear(a1);
      nmod_poly_clear(a2);
      nmod_poly_clear(b);
      nmod_poly_clear(c);
   }

   printf("PASS\n");
   return 0;
}
