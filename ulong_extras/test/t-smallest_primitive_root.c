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

    Copyright (C) 2013 Marcin Bodych

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include <math.h>

int main(void)
{
   int i,j, result;
   ulong k;
   mp_limb_t d, r, bits;
   mpz_t d_m;
   n_factor_t fac;
   mp_limb_signed_t exp1;
   mp_limb_t phi;
   flint_rand_t state;
   flint_randinit(state);

   printf("smallest_primitive_root....");
   fflush(stdout);
   for (j = 0; j < 10000 * flint_test_multiplier(); j++) 
   {
	
      mpz_init(d_m);
	  /* generating cases */
	  if( j == 0)
		 d = 2;
	  else if( j == 1)
	     d = 4;
	  else
	  {
		  /* less bigger powers */
		 k = 13 - (int)exp(log((double)(j%609))/2.5);
		 
		 do
		 {
			bits = n_randint(state, FLINT_D_BITS / k);
            d = n_randtest_bits(state, bits);
            mpz_set_ui(d_m, d);
            mpz_nextprime(d_m, d_m);
            d = mpz_get_ui( d_m );
		 } while (d < 3);
         d = n_pow(d, k ) * (1 + (j % 2));
         
	  }
	  
      r = n_smallest_primitive_root(d);
      result = 1;
      phi = n_euler_phi(d);
	  if(n_powmod(r, phi, d) != 1) 
	  {
         result = 0;
	  }
	  if(result)
	  {
		  n_factor_init(&fac);
		  n_factor(&fac, phi, 1);
		  
		  for (i = 0; i < fac.num; i++)
		  {
			 exp1 = phi / fac.p[i];
			 
			 if(n_powmod(r, exp1, d) == 1) 
			 {
				result = 0;
			 }
		  }
      }  
      if (!result)
      {
         printf("FAIL:\n");
         printf("d = %lu, r = %lu \n", d, r); 
         abort();
      }

      mpz_clear(d_m);
   }

   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
