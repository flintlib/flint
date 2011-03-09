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
   Copyright (C) 2010 Daniel Woodhouse

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int main(void)
{
   fmpz_t res;
   fmpz_t r1;
   fmpz_t m1;
   fmpz_t mm2;
   fmpz_t c;
   
   ulong r2 = 55947;
   ulong m2 = 104551;
   
   double pre = n_precompute_inverse(m2);

   printf("CRT....");
   fflush(stdout);

   fmpz_init(res);
   fmpz_init(r1);
   fmpz_init(m1);
   fmpz_init(mm2);
   fmpz_init(c);

   fmpz_set_ui(r1, 55953);
   fmpz_set_ui(m1, 104549);
   fmpz_set_ui(mm2, 104551);
     
   fmpz_invmod(c, m1, mm2);
   
   fmpz_CRT_ui_precomp(res, r1, m1, r2, m2, fmpz_get_ui(c), pre);
   
   if(fmpz_get_ui(res) != 369600)   
   {
      printf("FAIL");
      fmpz_print(res); printf("\n");
      abort();
   }

#if FLINT_BITS == 64
   fmpz_set_ui(r1, 0);
   fmpz_set_ui(m1, 18446744073709551557UL);
   fmpz_set_ui(mm2, 18446744073709551533UL);
   r2 = 0;
   m2 = 18446744073709551533UL;
    
   fmpz_set_ui(c, 14603672391686728297UL);
   pre = n_precompute_inverse(m2);

   fmpz_CRT_ui2_precomp(res, r1, m1, r2, m2, fmpz_get_ui(c), pre);
   if(fmpz_get_ui(res) != 0)
   {
      printf("FAIL");
      fmpz_print(res); printf("\n");
      abort();
   }
#endif

   fmpz_clear(r1);
   fmpz_clear(m1);
   fmpz_clear(mm2);
   fmpz_clear(res);
   fmpz_clear(c);

   printf("PASS\n");
   
   return 0;
}
