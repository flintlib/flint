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
   int i, j, k;
   mp_limb_t small_factor, p, pinv, pmod, x;
   fmpz_t n;
   fmpz_t temp, temp2, B, B_2;
   qs_t qs_inf;
   fmpz_init(n);
   fmpz_init(temp);
   fmpz_init(temp2);
   fmpz_init(B_2);
   FLINT_TEST_INIT(state);

   flint_printf("compute_poly_data....");
   fflush(stdout);

   for (i = 0; i < 1000; i++)
   {
       fmpz_randtest_unsigned(n, state, 130);

       if (fmpz_is_zero(n) || fmpz_is_one(n)) continue;

       qsieve_init(qs_inf, n);
       small_factor = qsieve_knuth_schroeppel(qs_inf);

       if (small_factor) continue;

       fmpz_mul_ui(qs_inf->kn, qs_inf->n, qs_inf->k); /* haven't calculated earlier */
       small_factor = qsieve_primes_init(qs_inf);

       if (small_factor) continue;

       qsieve_poly_init(qs_inf);

       small_factor = qsieve_compute_A(qs_inf);

       if (small_factor) continue;


       /* check if prime factor of 'A0', are it's factor */

       for (j = 0; j < qs_inf->s; j++)
       {
           if (fmpz_fdiv_ui(qs_inf->A0, qs_inf->factor_base[qs_inf->A_ind[j]].p) != 0)
           {
               abort();
           }
       }

       qsieve_compute_pre_data(qs_inf);

       /* check if precomputed values for 'A0' are correct */

       for (j = 0; j < qs_inf->s; j++)
       {
           p = qs_inf->factor_base[qs_inf->A_ind[j]].p;
           pinv = n_preinvert_limb(p);
           fmpz_divexact_ui(temp, qs_inf->A0, p);

           if (fmpz_cmp(temp, qs_inf->A0_divp[j]) != 0)
           {
               abort();
           }

           pmod = fmpz_fdiv_ui(qs_inf->A0_divp[j], p);
           x = n_invmod(pmod, p);
           x = n_mulmod2_preinv(x, qs_inf->sqrts[qs_inf->A_ind[j]], p, pinv);

           if (qs_inf->B0_terms[j] != x)
           {
               abort();
           }

       }

       /* check if all the 'B' values for particular hypercube 'A' satisfy,
          $B^2 \equiv kn \pmod{A}$
       */
       for (j = 0; j < qs_inf->num_q0; j++)
       {
           qs_inf->q0 = qs_inf->q0_values[j];
           fmpz_mul_ui(qs_inf->A, qs_inf->A0, qs_inf->q0);
           fmpz_mod(temp2, qs_inf->kn, qs_inf->A);
           qsieve_init_poly_first(qs_inf);
           qsieve_init_poly_next(qs_inf);

           for (k = 0; k < (1 << qs_inf->s); k++)
           {
               fmpz_pow_ui(B_2, qs_inf->B[k], 2);
               fmpz_mod(temp, B_2, qs_inf->A);

               if (fmpz_cmp(temp, temp2) != 0)
               {
                   abort();
               }

           }
       }

       qsieve_clear(qs_inf);
   }

   fmpz_clear(n);
   fmpz_clear(temp);
   fmpz_clear(temp2);
   fmpz_clear(B_2);

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}

