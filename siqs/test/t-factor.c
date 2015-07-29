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

    Copyright (C) 2015 Nitin Kumar

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

int main(void)
{
   int i, j = 0;
   qs_t qs_inf;
   mp_limb_t small_factor, a, b;
   fmpz_t n, x, y;
   fmpz_init(x);
   fmpz_init(y);
   fmpz_factor_t factors;
   FLINT_TEST_INIT(state);

   flint_printf("factor....");
   fflush(stdout);



   for (i = 0; i < 1; i++) /* Test random n */
   {
      b = a = n_randprime(state, 20, 1);

      while (b == a)
        b = n_randprime(state, 20, 1);

      fmpz_init_set_ui(n, a);
      fmpz_mul_ui(n, n, b);

      fmpz_set_str(x, "282174488599599500573849980909", 10);

      fmpz_set_str(y, "671998030559713968361666935769", 10);

      fmpz_mul(n, x, y);
      //flint_printf("\n n = %wu * %wu\n", a, b);
     // fmpz_set_str(n, "673899295409", 10);
      //fmpz_print(n);
      flint_printf("\n");

     // if (i == 0) continue;

     // fmpz_factor_init(factors);

      qsieve_factor(n, factors);

     // qsieve_process_partial(qs_inf);

     // fmpz_factor_clear(factors);

      //break;
   }

   fmpz_clear(n);
   fmpz_clear(x);
   fmpz_clear(y);

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}
