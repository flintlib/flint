/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qsieve.h"

void qsieve_square_root(fmpz_t X, fmpz_t Y, qs_t qs_inf,
   uint64_t * nullrows, slong ncols, slong l, fmpz_t N)
{
   slong position, i, j;
   slong * relation = qs_inf->relation;
   prime_t * factor_base = qs_inf->factor_base;
   slong * prime_count = qs_inf->prime_count;
   slong num_primes = qs_inf->num_primes;
   fmpz * Y_arr = qs_inf->Y_arr;
   fmpz_t pow;

   fmpz_init(pow);

   memset(prime_count, 0, num_primes*sizeof(slong));

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
#if QS_DEBUG
      if (prime_count[i] & 1)
      {
         flint_printf("Not a square in square root, %ld^%ld\n", factor_base[i].p, prime_count[i]);
         flint_printf("This is prime %ld of %ld factor base primes and %ld ks primes\n", i, qs_inf->num_primes, qs_inf->ks_primes);
      }
#endif
      if (prime_count[i])
      {
         fmpz_set_si(pow, factor_base[i].p);
         fmpz_powm_ui(pow, pow, prime_count[i]/2, N);
         fmpz_mul(X, X, pow);
      }

      if (i % 10 == 0 || i == num_primes - 1) fmpz_mod(X, X, N);
   }

   fmpz_clear(pow);
}

