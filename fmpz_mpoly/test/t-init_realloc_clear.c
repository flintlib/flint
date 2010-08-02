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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int main(void)
{
   int i, result;
   printf("init/init2/realloc/clear....");
   fflush(stdout);
   
   for (i = 0; i < 10000; i++) 
   {
      fmpz_mpoly_t a;

      fmpz_mpoly_init2(a, n_randint(100), n_randint(10), n_randint(100));
      fmpz_mpoly_clear(a);      
   }

   for (i = 0; i < 10000; i++) 
   {
      fmpz_mpoly_t a;

      fmpz_mpoly_init2(a, n_randint(100), n_randint(10), n_randint(100));
      fmpz_mpoly_realloc(a, n_randint(100));
      fmpz_mpoly_realloc(a, n_randint(100));
      fmpz_mpoly_clear(a);      
   }
      
   _fmpz_cleanup();
   printf("PASS\n");
   return 0;
}
