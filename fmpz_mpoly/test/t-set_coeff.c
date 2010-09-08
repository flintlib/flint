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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"
#include <time.h>

int main(void)
{
   int result;
   printf("set_coeff....");
   fflush(stdout);

   fmpz_mpoly_t poly1;
   fmpz_mpoly_t poly2;

   fmpz_t temp;
   fmpz_init(temp);

   fmpz_mpoly_init2(poly1, 10, 3, 21);
   fmpz_mpoly_init2(poly2, 10, 3, 21);

   fmpz_set_ui(poly1->coeffs + 0, 1);
   fmpz_set_ui(poly1->coeffs + 1, 1);
   fmpz_set_ui(poly1->coeffs + 2, 1);

   fmpz_set_ui(poly1->exps + 0, 1);
   fmpz_set_ui(poly1->exps + 1, (ulong)1 << 21);
   fmpz_set_ui(poly1->exps + 2, (ulong)1 << 42);

   poly1->length = 3;

   //now set poly2
   fmpz_set_ui(temp, 1);
   fmpz_mpoly_set_coeff_fmpz(poly2, (ulong)1 <<42, temp);
   fmpz_mpoly_set_coeff_fmpz(poly2, (ulong)1 <<21, temp);
   fmpz_mpoly_set_coeff_fmpz(poly2, 1, temp);

   if(fmpz_mpoly_equal(poly1, poly2) == 0){
      printf("FAIL");
   }

 

   fmpz_mpoly_set_coeff_fmpz(poly2, 3, temp);
   fmpz_set_ui(temp, 0);
   fmpz_mpoly_set_coeff_fmpz(poly2, 3, temp);
   
 
   if(fmpz_mpoly_equal(poly1, poly2) == 0){
      printf("FAIL");
   }
   
   printf("PASS\n");
   return 0;
}
