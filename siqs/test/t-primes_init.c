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

/******************************IGNORE*****************************************/

#include "C:\Users\measure\Documents\GitHub\flint2\siqs\qsieve.h"
#include "C:\Users\measure\Documents\GitHub\flint2\siqs\clear.c"
#include "C:\Users\measure\Documents\GitHub\flint2\siqs\init.c"
#include "C:\Users\measure\Documents\GitHub\flint2\siqs\knuth_schroeppel.c"
#include "C:\Users\measure\Documents\GitHub\flint2\siqs\primes_init.c"

/******************************************************************************/

int main(void)
{
   int i;
   slong j, k;
   mp_limb_t small_factor, pmod;
   qs_t qs_inf;
   fmpz_t n, x, y;
   fmpz_init(n);
   fmpz_init(x);
   fmpz_init(y);

   FLINT_TEST_INIT(state);

   flint_printf("primes_init....");
   fflush(stdout);

   for (i = 0; i < 100000; i++)
   {

       fmpz_randtest_unsigned(n, state, 130);

       if (fmpz_is_zero(n) || fmpz_is_one(n)) continue;

       qsieve_init(qs_inf, n);
       small_factor = qsieve_knuth_schroeppel(qs_inf);

       if (small_factor) continue;

       fmpz_mul_ui(qs_inf->kn, n, qs_inf->k);                      /*haven't calculated earlier*/
       small_factor = qsieve_primes_init(qs_inf);

       if (small_factor)
       {
           /* check if factor, if returned by factor base function
              is actually a factor of n */

           if (fmpz_fdiv_ui(n, small_factor))
           {
               abort();                                       /*exception to add*/
           }
           else continue;
       }

       for (j = 2; j < qs_inf->num_primes; j++)
       {
           /* check if square root of kn modulo
               factor base prime p are correct   */

           fmpz_set_si(x, qs_inf->sqrts[j]);
           fmpz_mul(y, x, x);
           fmpz_mod_ui(x, y, qs_inf->factor_base[j].p);
           fmpz_mod_ui(y, qs_inf->kn, qs_inf->factor_base[j].p);

           if (fmpz_cmp(x, y) != 0)
           {
               abort();                                                    /*exception to add*/
           }

           /* check if inverse of factor base primes are correct*/

           fmpz_randtest_unsigned(x, state, FLINT_BITS);
           fmpz_mod_ui(y, x, qs_inf->factor_base[j].p);

           pmod = n_mod2_preinv(fmpz_get_ui(x),
                        qs_inf->factor_base[j].p, qs_inf->factor_base[j].pinv);

           if (fmpz_get_ui(y) != pmod)
           {
               abort();                                                    /*exception to add*/
           }
       }

       /* add next 50 primes to factor base */

       k = qs_inf->num_primes;

       small_factor = qsieve_primes_increment(qs_inf, 50);

       if (small_factor)
       {
           if (fmpz_fdiv_ui(qs_inf->n, small_factor))
           {
               abort();                                                   /* exception to add */
           }
           else continue;
       }


       for (j = k; j < qs_inf->num_primes; j++)
       {

           fmpz_set_si(x, qs_inf->sqrts[j]);
           fmpz_mul(y, x, x);
           fmpz_mod_ui(x, y, qs_inf->factor_base[j].p);
           fmpz_mod_ui(y, qs_inf->kn, qs_inf->factor_base[j].p);

           if (fmpz_cmp(x, y) != 0)
           {
               abort();
           }


           fmpz_randtest_unsigned(x, state, FLINT_BITS);
           fmpz_mod_ui(y, x, qs_inf->factor_base[j].p);

           pmod = n_mod2_preinv(fmpz_get_ui(x),
                        qs_inf->factor_base[j].p, qs_inf->factor_base[j].pinv);

           if (fmpz_get_ui(y) != pmod)
           {
               abort();
           }
       }

       /* add next 30 primes to factor base */

       k = qs_inf->num_primes;

       small_factor = qsieve_primes_increment(qs_inf, 30);

       if (small_factor)
       {
           if (fmpz_fdiv_ui(qs_inf->n, small_factor))
           {
               abort();                                         /* exception to add */
           }
           else continue;
       }


       for (j = k; j < qs_inf->num_primes; j++)
       {

           fmpz_set_si(x, qs_inf->sqrts[j]);
           fmpz_mul(y, x, x);
           fmpz_mod_ui(x, y, qs_inf->factor_base[j].p);
           fmpz_mod_ui(y, qs_inf->kn, qs_inf->factor_base[j].p);

           if (fmpz_cmp(x, y) != 0)
           {
               abort();
           }


           fmpz_randtest_unsigned(x, state, FLINT_BITS);
           fmpz_mod_ui(y, x, qs_inf->factor_base[j].p);

           pmod = n_mod2_preinv(fmpz_get_ui(x),
                        qs_inf->factor_base[j].p, qs_inf->factor_base[j].pinv);

           if (fmpz_get_ui(y) != pmod)
           {
               abort();
           }
       }


       /* add next 10 primes to factor base*/

       k = qs_inf->num_primes;

       small_factor = qsieve_primes_increment(qs_inf, 10);

       if (small_factor)
       {
           if (fmpz_fdiv_ui(qs_inf->n, small_factor))
           {
               abort();                                         /* exception to add */
           }
           else continue;
       }


       for (j = k; j < qs_inf->num_primes; j++)
       {

           fmpz_set_si(x, qs_inf->sqrts[j]);
           fmpz_mul(y, x, x);
           fmpz_mod_ui(x, y, qs_inf->factor_base[j].p);
           fmpz_mod_ui(y, qs_inf->kn, qs_inf->factor_base[j].p);

           if (fmpz_cmp(x, y) != 0)
           {
               abort();
           }


           fmpz_randtest_unsigned(x, state, FLINT_BITS);
           fmpz_mod_ui(y, x, qs_inf->factor_base[j].p);

           pmod = n_mod2_preinv(fmpz_get_ui(x),
                        qs_inf->factor_base[j].p, qs_inf->factor_base[j].pinv);

           if (fmpz_get_ui(y) != pmod)
           {
               abort();
           }
       }



       qsieve_clear(qs_inf);
   }

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}

