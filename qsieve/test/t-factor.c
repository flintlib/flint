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
    Copyright (C) 2016 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

void randprime(fmpz_t p, flint_rand_t state, slong bits)
{
    fmpz_randbits(p, state, bits);
 
    if (fmpz_sgn(p) < 0)
       fmpz_neg(p, p);

    if (fmpz_is_even(p))
       fmpz_add_ui(p, p, 1);
 
    while (!fmpz_is_probabprime(p))
       fmpz_add_ui(p, p, 2);
}

int main(void)
{
   slong i;
   fmpz_t n, x, y, z;
   fmpz_factor_t factors;
   FLINT_TEST_INIT(state);

   fmpz_init(x);
   fmpz_init(y);
   fmpz_init(z);
   fmpz_init(n);

   flint_printf("factor....");
   fflush(stdout);

   /* Test n with large prime factor */
   {
      fmpz_set_str(n, "12387192837918273918723981291837121933111751252512531193171", 10);
    
      fmpz_factor_init(factors);

      qsieve_factor(factors, n);

      if (factors->num < 5)
      {
         flint_printf("FAIL:\n");
         flint_printf("%ld factors found\n", factors->num);
         abort();
      }

      fmpz_factor_clear(factors);
   }

   for (i = 0; i < 30; i++) /* Test random n, two factors */
   {
      slong bits = 40;

      randprime(x, state, bits);
      do {
         randprime(y, state, bits);
      } while (fmpz_equal(x, y));
      
      fmpz_mul(n, x, y);

      fmpz_factor_init(factors);

      qsieve_factor(factors, n);

      if (factors->num < 2)
      {
         flint_printf("FAIL:\n");
         flint_printf("%ld factors found\n", factors->num);
         abort();
      }

      fmpz_factor_clear(factors);
   }

   for (i = 0; i < 30; i++) /* Test random n, three factors */
   {
      randprime(x, state, 40);
      do {
         randprime(y, state, 40);
      } while (fmpz_equal(x, y));
      do {
         randprime(z, state, 40);
      } while (fmpz_equal(x, z) || fmpz_equal(y, z));

      fmpz_mul(n, x, y);
      fmpz_mul(n, n, z);

      fmpz_factor_init(factors);

      qsieve_factor(factors, n);

      if (factors->num < 3)
      {
         flint_printf("FAIL:\n");
         flint_printf("%ld factors found\n", factors->num);
         abort();
      }

      fmpz_factor_clear(factors);
   }

   for (i = 0; i < 30; i++) /* Test random n, small factors */
   {
      randprime(x, state, 10);
      do {
         randprime(y, state, 10);
      } while (fmpz_equal(x, y));
      randprime(z, state, 40);

      fmpz_mul(n, x, y);
      fmpz_mul(n, n, z);

      fmpz_factor_init(factors);

      qsieve_factor(factors, n);

      if (factors->num < 3)
      {
         flint_printf("FAIL:\n");
         flint_printf("%ld factors found\n", factors->num);
         abort();
      }

      fmpz_factor_clear(factors);
   }

   fmpz_clear(n);
   fmpz_clear(x);
   fmpz_clear(y);
   fmpz_clear(z);

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}
