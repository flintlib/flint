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
   int i;
   mp_limb_t small_factor;
   fmpz_t n;
   qs_t qs_inf;
   fmpz_init(n);
   FLINT_TEST_INIT(state);

   flint_printf("compute_poly_data....");
   fflush(stdout);

   for (i = 0; i < 100; i++)
   {
       fmpz_randtest_unsigned(n, state, 130);

       if (fmpz_is_zero(n) || fmpz_is_one(n)) continue;

       qsieve_init(qs_inf, n);
       small_factor = qsieve_knuth_schroeppel(qs_inf);
       if (small_factor) continue;
       fmpz_mul_ui(qs_inf->kn, qs_inf->n, qs_inf->k);
       small_factor = qsieve_primes_init(qs_inf);
       if (small_factor) continue;
       flint_printf("Number to be factored:\n");
       fmpz_print(qs_inf->kn);
       flint_printf("\n");

       flint_printf("Optimal value of A: %wd \n", qs_inf->target_A);

       qsieve_compute_A0(qs_inf);

       flint_printf("Approximation to target value: %wd \n", qs_inf->A0);
   }

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}

