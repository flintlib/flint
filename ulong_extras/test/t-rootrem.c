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

    Copyright (C) 2015 Kushagra Singh

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
   
   flint_printf("rootrem....");
   fflush(stdout);

   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t a, b, c, d, i, j;
      mpz_t e, f, g, h;
      int res;

      mpz_init(e);
      mpz_init(f);
      mpz_init(g);
      mpz_init(h);
      
      c = n_randint(state, 0);    /*number */
      flint_mpz_set_ui(g, c);


      d = n_randint(state, 0);   /*root */
      flint_mpz_set_ui(h, d);

      res = n_rootrem(&a, &b, c, d);

      mpz_rootrem(e, f, g, flint_mpz_get_ui(h));
      
      i = flint_mpz_get_ui(e);
      j = flint_mpz_get_ui(f);

      result = ((a == i) && (b == j));

      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("Passed Parameters : n = %wu root = %wu", c, d);
         flint_printf("Answer generated : base = %wu remainder = %wu", a, b);
         flint_printf("Expected answer : base = %wu remainder = %wu", i, j);
         abort();
      }


      mpz_clear(e);
      mpz_clear(f);
      mpz_clear(g);
      mpz_clear(h);
   }


   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}