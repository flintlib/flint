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

    Copyright (C) 2006, 2011 William Hart

******************************************************************************/

#undef ulong /* avoid clash with stdlib */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long 

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"
#include "fmpz.h"

void qsieve_ll_square_root(fmpz_t X, fmpz_t Y, qs_t qs_inf, 
   uint64_t * nullrows, long ncols, long l, fmpz_t N)
{
   long position, i, j;
   long * relation = qs_inf->relation;
   prime_t * factor_base = qs_inf->factor_base;
   long * prime_count = qs_inf->prime_count;
   long num_primes = qs_inf->num_primes;
   fmpz * Y_arr = qs_inf->Y_arr; 
   fmpz_t pow;

   fmpz_init(pow);
      
   memset(prime_count, 0, num_primes*sizeof(long));
      
   fmpz_one(X);
   fmpz_one(Y);
   
   for (i = 0; i < ncols; i++)
   {
      if (get_null_entry(nullrows, i, l)) 
      {
         position = qs_inf->matrix[i].orig*2*qs_inf->max_factors;
         
         for (j = 0; j < relation[position]; j++)
         {
            prime_count[relation[position+2*j+1]] +=
               (relation[position+2*j+2]);
         }
         
         fmpz_mul(Y, Y, Y_arr + qs_inf->matrix[i].orig);
         if (i % 10 == 0) fmpz_mod(Y, Y, N);
      }
   }

   for (i = 0; i < num_primes; i++)
   {
      if (prime_count[i]) 
      {
         fmpz_set_ui(pow, factor_base[i].p);
         fmpz_powm_ui(pow, pow, prime_count[i]/2, N);
         fmpz_mul(X, X, pow);
      } 

      if (i%10 == 0 || i == num_primes - 1) fmpz_mod(X, X, N);
   }

#if QS_DEBUG
   for (i = 0; i < num_primes; i++)
      if ((prime_count[i] %2) != 0) printf("Error %ld, %ld, %ld\n", l, i, prime_count[i]);
#endif

   fmpz_clear(pow);
}

