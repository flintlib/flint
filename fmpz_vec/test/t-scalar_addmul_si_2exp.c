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

   Copyright (C) 2009 William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int main(void)
{
   int result;
   printf("scalar_addmul_si_2exp....");
   fflush(stdout);
   
   fmpz_poly_randinit();
   
   // compare with alternative method of computation
   for (ulong i = 0; i < 10000UL; i++) 
   {
      fmpz_poly_t a, b, c;
      ulong length = n_randint(100);
      long x;
	  mp_bitcnt_t exp;

      fmpz_poly_init2(a, length);
      fmpz_poly_init2(b, length);
      fmpz_poly_init2(c, length);
	  
      fmpz_poly_randtest(a, length, n_randint(200));
      fmpz_poly_randtest(b, length, n_randint(200));
      fmpz_poly_set(c, b);
	  x = n_randbits(n_randint(FLINT_BITS - 1));
      if (n_randint(2)) x = -x;
      exp = n_randint(200);

      _fmpz_vec_scalar_addmul_si_2exp(b->coeffs, a->coeffs, length, x, exp);
      _fmpz_vec_scalar_addmul_si_2exp(c->coeffs, a->coeffs, length, x, exp);
      
      result = (fmpz_poly_equal(b, c));
      if (!result)
      {
         printf("Error:\n");
         printf("x = %ld, exp = %lu\n", x, exp);
		 fmpz_poly_print(b); printf("\n\n");
         fmpz_poly_print(c); printf("\n\n");
         abort();
      }

      fmpz_poly_clear(a);
      fmpz_poly_clear(b);
      fmpz_poly_clear(c);
   }

   fmpz_poly_randclear();
      
   _fmpz_cleanup();
   printf("PASS\n");
   return 0;
}
