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

    Copyright (C) 2009, 2011 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

/******************************IGNORE*****************************************/

#include "C:\Users\measure\Documents\GitHub\flint2\siqs\qsieve.h"
#include "C:\Users\measure\Documents\GitHub\flint2\siqs\clear.c"
#include "C:\Users\measure\Documents\GitHub\flint2\siqs\init.c"
#include "C:\Users\measure\Documents\GitHub\flint2\siqs\knuth_schroeppel.c"

/******************************************************************************/

int main(void)
{
   int i;
   mp_limb_t small_factor;
   FLINT_TEST_INIT(state);

   flint_printf("knuth_schroeppel....");
   fflush(stdout);



   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random n */
   {
      fmpz_t n;
      fmpz_init(n);
      qs_t qs_inf;

      fmpz_randtest_unsigned(n, state, 130);

      if (fmpz_is_zero(n) || fmpz_is_one(n)) continue;

      qsieve_init(qs_inf, n);
      small_factor = qsieve_knuth_schroeppel(qs_inf);

      if (small_factor)
      {
          if (fmpz_fdiv_ui(qs_inf->n, small_factor))
          {
              flint_printf("%wd is not a factor of ", small_factor);
              fmpz_print(qs_inf->n);
              flint_printf("\n");
              abort();
          }
      }

      qsieve_clear(qs_inf);
   }

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}
