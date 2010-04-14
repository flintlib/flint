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

   Copyright (C) 2009, 2010 William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int main(void)
{
   int result;
   printf("scalar_mul_2exp....");
   fflush(stdout);
   
   _fmpz_vec_randinit();
   
   // check aliasing of a and b
   for (ulong i = 0; i < 10000UL; i++) 
   {
      fmpz * a, * b, * c;
      ulong length = n_randint(100);
      ulong exp = n_randint(200);
      
      a = _fmpz_vec_init(length);
      b = _fmpz_vec_init(length);
      _fmpz_vec_randtest(a, length, n_randint(200));
      
      _fmpz_vec_scalar_mul_2exp(b, a, length, exp);
      _fmpz_vec_scalar_mul_2exp(a, a, length, exp);
      
      result = (_fmpz_vec_equal(a, b, length));
      if (!result)
      {
         printf("Error:\n");
         printf("exp = %lu\n", exp);
		 _fmpz_vec_print(a, length); printf("\n\n");
         _fmpz_vec_print(b, length); printf("\n\n");
         abort();
      }

      _fmpz_vec_clear(a, length);
      _fmpz_vec_clear(b, length);
   }

   // check aliasing of (a*2^e1)*2^e2 equals a*2^(e1+e2)
   for (ulong i = 0; i < 10000UL; i++) 
   {
      fmpz * a, * b, * c;
      ulong length = n_randint(100);
      ulong e1 = n_randint(200);
      ulong e2 = n_randint(200);
      
      a = _fmpz_vec_init(length);
      b = _fmpz_vec_init(length);
      _fmpz_vec_randtest(a, length, n_randint(200));
      
      _fmpz_vec_scalar_mul_2exp(b, a, length, e1);
      _fmpz_vec_scalar_mul_2exp(b, b, length, e2);
      _fmpz_vec_scalar_mul_2exp(a, a, length, e1 + e2);
      
      result = (_fmpz_vec_equal(a, b, length));
      if (!result)
      {
         printf("Error:\n");
         printf("e1 = %lu, e2 = %lu\n", e1, e2);
		 _fmpz_vec_print(a, length); printf("\n\n");
         _fmpz_vec_print(b, length); printf("\n\n");
		 abort();
      }

      _fmpz_vec_clear(a, length);
      _fmpz_vec_clear(b, length);
   }

   _fmpz_vec_randclear();
      
   _fmpz_cleanup();
   printf("PASS\n");
   return 0;
}
