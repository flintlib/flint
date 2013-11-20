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
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int i, result;
   FLINT_TEST_INIT(state);
   

   flint_printf("add_ssaaaa....");
   fflush(stdout);

   for (i = 0; i < 1000000; i++)
   {
      mp_limb_t sh1, sl1, sh2, sl2, ah1, al1, ah2, al2;

      ah1 = n_randtest(state);
      al1 = n_randtest(state);
      ah2 = n_randtest(state);
      al2 = n_randtest(state);
      
      add_ssaaaa(sh1, sl1, ah1, al1, ah2, al2);
      
      sl2 = al1 + al2;
      sh2 = (sl1 < al1);
      sh2 += ah1;
      sh2 += ah2;

      result = ((sh2 == sh1) && (sl2 == sl1));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("ah1 = %wu, al1 = %wu, ah2 = %wu, al1 = %wu\n", ah1, al1, ah2, al1); 
         flint_printf("sh2 = %wu, sh1 = %wu, sl2 = %wu, sl1 = %wu\n", sh2, sh1, sl2, sl1);
         abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
