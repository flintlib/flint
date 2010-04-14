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
#include "fmpz_vec.h"
#include "ulong_extras.h"

int main(void)
{
   int result;
   printf("init/clear....");
   fflush(stdout);
   
   // compare with alternative method of computation
   for (ulong i = 0; i < 10000UL; i++) 
   {
      fmpz * a;
      ulong j;
	  
	  ulong length = n_randint(100) + 1;
      
	  a = _fmpz_vec_init(length);
      for (j = 0; j < length; j++)
	     fmpz_zero(a + j);

	  _fmpz_vec_clear(a, length);
   }

   _fmpz_cleanup();
   printf("PASS\n");
   return 0;
}
