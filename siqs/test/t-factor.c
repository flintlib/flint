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

    /*  fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

    /*  fmpz_set_str(x, "282174488599599500573849980909", 10);

      fmpz_set_str(y, "280829369862134719390036617067", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

    /*  fmpz_set_str(x, "513821217024129243948411056803", 10);

      fmpz_set_str(y, "521419622856657689423872613771", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

   /*   fmpz_set_str(x, "341427877364219557396646723561", 10);

      fmpz_set_str(y, "359916012598740083996400089999", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

   /*   fmpz_set_str(x, "564819669946735512444543556507", 10);

      fmpz_set_str(y, "590872612825179551336102196593", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

   /*   fmpz_set_str(x, "341427877364219557396646723561", 10);

      fmpz_set_str(y, "359916012598740083996400089999", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

  /*    fmpz_set_str(x, "666666777777788888888999999999", 10);

      fmpz_set_str(y, "671998030559713968361666935769", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

      fmpz_set_str(x, "1291713871775387910289126481863", 10);

      fmpz_set_str(y, "1317313773719779779173773137131", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors);

   /*   fmpz_set_str(x, "3111111111111111111111111111113", 10);

      fmpz_set_str(y, "3113113113113113113113113113113", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

   /*   fmpz_set_str(x, "1011313133831810181383313131101", 10);

      fmpz_set_str(y, "1124495649695212998778414534811", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

   /*   fmpz_set_str(x, "909090909090909090909090909091", 10);

      fmpz_set_str(y, "811451682377384625400019885321", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

   /*   fmpz_set_str(x, "1000000000000000000000000000057", 10);

      fmpz_set_str(y, "1000000000000066600000000000001", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */

   /*   fmpz_set_str(x, "1717171717171717171717171717171", 10);

      fmpz_set_str(y, "2771826449193354891007108898387", 10);

      fmpz_mul(n, x, y);

      fmpz_print(n);

      flint_printf("\n");

      qsieve_factor(n, factors); */
   }

   fmpz_clear(n);
   fmpz_clear(x);
   fmpz_clear(y);

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}
