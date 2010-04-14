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
   printf("sub....");
   fflush(stdout);
   
   _fmpz_vec_randinit();
   
   // check aliasing of a and c
   for (ulong i = 0; i < 10000UL; i++) 
   {
      fmpz * a, * b, * c;
      ulong length = n_randint(100);
      
      a = _fmpz_vec_init(length);
      b = _fmpz_vec_init(length);
      c = _fmpz_vec_init(length);
      _fmpz_vec_randtest(a, length, n_randint(200));
      _fmpz_vec_randtest(b, length, n_randint(200));
      
      _fmpz_vec_sub(c, a, b, length);
      _fmpz_vec_sub(a, a, b, length);
      
      result = (_fmpz_vec_equal(a, c, length));
      if (!result)
      {
         printf("Error:\n");
         _fmpz_vec_print(a, length); printf("\n\n");
         _fmpz_vec_print(c, length); printf("\n\n");
         abort();
      }

      _fmpz_vec_clear(a, length);
      _fmpz_vec_clear(b, length);
      _fmpz_vec_clear(c, length);
   }

   // check aliasing of b and c
   for (ulong i = 0; i < 10000UL; i++) 
   {
      fmpz * a, * b, * c;
      ulong length = n_randint(100);
      
      a = _fmpz_vec_init(length);
      b = _fmpz_vec_init(length);
      c = _fmpz_vec_init(length);
      _fmpz_vec_randtest(a, length, n_randint(200));
      _fmpz_vec_randtest(b, length, n_randint(200));
      
      _fmpz_vec_sub(c, a, b, length);
      _fmpz_vec_sub(b, a, b, length);
      
      result = (_fmpz_vec_equal(b, c, length));
      if (!result)
      {
         printf("Error:\n");
         _fmpz_vec_print(b, length); printf("\n\n");
         _fmpz_vec_print(c, length); printf("\n\n");
         abort();
      }

      _fmpz_vec_clear(a, length);
      _fmpz_vec_clear(b, length);
      _fmpz_vec_clear(c, length);
   }

   // check a + b - b = a
   for (ulong i = 0; i < 10000UL; i++) 
   {
      fmpz * a, * b, * c, * d;
      ulong length = n_randint(100);
      
      a = _fmpz_vec_init(length);
      b = _fmpz_vec_init(length);
      c = _fmpz_vec_init(length);
      d = _fmpz_vec_init(length);
      _fmpz_vec_randtest(a, length, n_randint(200));
      _fmpz_vec_randtest(b, length, n_randint(200));
      
      _fmpz_vec_add(c, a, b, length);
      _fmpz_vec_sub(d, c, b, length);
      
      result = (_fmpz_vec_equal(d, a, length));
      if (!result)
      {
         printf("Error:\n");
         _fmpz_vec_print(a, length); printf("\n\n");
         _fmpz_vec_print(d, length); printf("\n\n");
         abort();
      }

      _fmpz_vec_clear(a, length);
      _fmpz_vec_clear(b, length);
      _fmpz_vec_clear(c, length);
      _fmpz_vec_clear(d, length);
   }


   _fmpz_vec_randclear();
      
   _fmpz_cleanup();
   printf("PASS\n");
   return 0;
}
