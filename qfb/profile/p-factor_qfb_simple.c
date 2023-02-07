/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright 2012 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "fmpz.h"
#include "qfb.h"

int main(void)
{
   fmpz_t n, p;
   slong iters, i, j;
   qfb_t pow, oldpow, twopow;
   ulong pr, oldpr, nmodpr;
   int done = 1;

   fmpz_init(n);
   fmpz_init(p);

   qfb_init(pow);
   qfb_init(oldpow);
   qfb_init(twopow);
   
   printf("Enter number to be factored: "); fflush(stdout);
   if (!fmpz_read(n))
   {
      printf("Read failed\n");
      abort();
   }
   fmpz_neg(n, n);

   printf("Enter a number of iterations: "); fflush(stdout);
   if (!scanf("%ld", &iters))
   {
      printf("Read failed\n");
      abort();
   }
    
   /* find prime such that n is a square mod p (or p divides n) */
   if (fmpz_is_even(n))
   {
      printf("Factor: 2\n");
      return 0;
   }
   
   pr = 2;
   
   if (fmpz_fdiv_ui(n, 4) == 3)
   {
      if (fmpz_fdiv_ui(n, 3) == 0)
      {
         printf("Factor: 3\n");
         return 0;
      }
      fmpz_mul_ui(n, n, 3);
      pr = 3;
   }

   do
   {
      pr = n_nextprime(pr, 0);
      nmodpr = fmpz_fdiv_ui(n, pr);
      
      if (nmodpr == 0) /* pr is a factor */
      {
         printf("Factor: %lu\n", pr);
         return 0;
      }
   } while (n_jacobi(nmodpr, pr) < 0);

   fmpz_set_ui(p, pr);

   /* find prime form of discriminant n */
   qfb_prime_form(pow, n, p);
   
   /* raise to various powers of small primes */
   qfb_pow_ui(pow, pow, n, 59049); /* 3^10 */
   qfb_pow_ui(twopow, pow, n, 4096);
   if (qfb_is_principal_form(twopow, n))
      goto done;

   qfb_pow_ui(pow, pow, n, 390625); /* 5^8 */
   qfb_pow_ui(twopow, pow, n, 4096);
   if (qfb_is_principal_form(twopow, n))
      goto done;

   qfb_pow_ui(pow, pow, n, 117649); /* 7^6 */
   qfb_pow_ui(twopow, pow, n, 4096);
   if (qfb_is_principal_form(twopow, n))
      goto done;

   qfb_pow_ui(pow, pow, n, 14641); /* 11^4 */
   qfb_pow_ui(twopow, pow, n, 4096);
   if (qfb_is_principal_form(twopow, n))
      goto done;

   qfb_pow_ui(pow, pow, n, 169); /* 13^2 */
   qfb_pow_ui(twopow, pow, n, 4096);
   if (qfb_is_principal_form(twopow, n))
      goto done;

   pr = 13;
   for (i = 0; i < iters; )
   {
      j = FLINT_MIN(i + 1024, iters);
      qfb_set(oldpow, pow);
      oldpr = pr;
      for ( ; i < j; i++)
      {
         pr = n_nextprime(pr, 0);
         qfb_pow_ui(pow, pow, n, pr*pr);
      }
      qfb_pow_ui(twopow, pow, n, 4096);
      if (qfb_is_principal_form(twopow, n))
      {
         qfb_set(pow, oldpow);
         pr = oldpr;
         while (1)
         {
            pr = n_nextprime(pr, 0);
            qfb_pow_ui(pow, pow, n, pr*pr);
            qfb_pow_ui(twopow, pow, n, 4096);
            if (qfb_is_principal_form(twopow, n))
               goto done;
         }
      }
   }
   if (i == iters)
      done = 0;

done:
   if (done)
   {
      for (i = 0; i < 12; i++)
      {
         qfb_pow_ui(twopow, pow, n, 2);
         if (qfb_is_principal_form(twopow, n))
         {
            qfb_print(pow); printf("\n");
            return 0;
         }
         qfb_set(pow, twopow);
      }
   } else
      printf("Failed\n");

   qfb_clear(pow);
   qfb_clear(oldpow);
   qfb_clear(twopow);
   
   fmpz_clear(n);
   fmpz_clear(p);

   return 0;
}
