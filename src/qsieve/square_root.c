/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "longlong.h"
#include "fmpz.h"
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

   /*
      X = prod_i p_i^(e_i) mod N, where e_i = prime_count[i]/2.

      Writing the exponents in binary,

          prod_i p_i^(e_i) = prod_k Q_k^(2^k),

      where Q_k is the product of those p_i whose exponent has bit k set.
      Evaluating that by Horner from the top bit down costs one modular
      squaring per bit of the largest exponent, plus one word multiplication
      per prime, in place of one modular exponentiation per prime.  Each
      Q_k is accumulated in a single word and flushed into X only when the
      next multiplication would overflow.

      Note that factor_base[2].p is -1: that entry carries the sign of the
      polynomial value rather than a prime.  Its magnitude contributes
      nothing to the product, so the sign is accumulated as a parity and
      applied at the end.
   */
   {
      slong maxexp = 0, nbits, k;
      int neg = 0;

      for (i = 0; i < num_primes; i++)
      {
#if QS_DEBUG
         if (prime_count[i] & 1)
         {
            flint_printf("Not a square in square root, %ld^%ld\n", factor_base[i].p, prime_count[i]);
            flint_printf("This is prime %ld of %ld factor base primes and %ld ks primes\n", i, qs_inf->num_primes, qs_inf->ks_primes);
         }
#endif
         if (factor_base[i].p < 0 && ((prime_count[i] >> 1) & 1))
            neg ^= 1;

         /* entries with |p| <= 1 contribute nothing to the product, so
            they must not inflate the number of squarings either */
         if (FLINT_ABS(factor_base[i].p) > 1 && (prime_count[i] >> 1) > maxexp)
            maxexp = prime_count[i] >> 1;
      }

      nbits = maxexp ? (slong) (FLINT_BITS - flint_clz((ulong) maxexp)) : 0;

      for (k = nbits - 1; k >= 0; k--)
      {
         ulong acc = 1;

         fmpz_mul(X, X, X);
         fmpz_mod(X, X, N);

         for (i = 0; i < num_primes; i++)
         {
            if (((prime_count[i] >> 1) >> k) & 1)
            {
               slong sp = factor_base[i].p;
               ulong pr = (ulong) FLINT_ABS(sp);

               if (pr <= 1)   /* the sign entry, or an absent multiplier */
                  continue;

               if (acc > (~UWORD(0)) / pr)   /* would overflow a word */
               {
                  fmpz_mul_ui(X, X, acc);
                  fmpz_mod(X, X, N);
                  acc = 1;
               }

               acc *= pr;
            }
         }

         if (acc != 1)
         {
            fmpz_mul_ui(X, X, acc);
            fmpz_mod(X, X, N);
         }
      }

      if (neg)
      {
         fmpz_neg(X, X);
         fmpz_mod(X, X, N);
      }
   }
}
