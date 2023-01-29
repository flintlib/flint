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
   fmpz_t g, n, n0, p, L;
   slong iters, i, j, depth, jmax = 10;
   qfb_t pow, twopow;
   ulong pr, nmodpr, mult, n0mod4;
   qfb_hash_t * qhash;
   int done;

   fmpz_init(g);
   fmpz_init(n);
   fmpz_init(n0);
   fmpz_init(p);
   fmpz_init(L);

   qfb_init(pow);
   qfb_init(twopow);
   
   printf("Enter number to be factored: "); fflush(stdout);
   if (!fmpz_read(n0))
   {
      printf("Read failed\n");
      abort();
   }
   fmpz_neg(n0, n0);

   iters = 10;
    
   /* find prime such that n is a square mod p (or p divides n) */
   if (fmpz_is_even(n0))
   {
      printf("Factor: 2\n");
      return 0;
   }
   
   mult = 1;
   
   while (1) /* keep increasing iterations for each multiplier until done */
   {
      printf("iters = %ld, multipliers = %ld\n", iters, jmax);
   
      for (j = 0; j < jmax; j++) /* loop over jmax different multipliers */
      {

         done = 1;

         if (mult != 1 && fmpz_fdiv_ui(n0, mult) == 0)
         {
            printf("Factor: %ld\n", mult);
            return 0;
         }
      
         fmpz_mul_ui(n, n0, mult);
         if (fmpz_fdiv_ui(n, 4) == 3)
            fmpz_mul_2exp(n, n, 2);

         pr = 2;
         fmpz_abs(L, n);
         fmpz_root(L, L, 4);

         do
         {
            pr = n_nextprime(pr, 0);
            while (mult % pr == 0)
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

         depth = FLINT_BIT_COUNT(iters) + 1;
         qhash = qfb_hash_init(depth);

         pr = 13;
         for (i = 0; i < iters; i++)
         {
            pr = n_nextprime(pr, 0);
            qfb_pow_ui(pow, pow, n, pr*pr);
            qfb_pow_ui(twopow, pow, n, 4096);
            if (qfb_is_principal_form(twopow, n)) /* found factor */
               break;
            qfb_hash_insert(qhash, twopow, pow, i, depth);
         }

         if (i == iters) /* stage 2 */
         {
            ulong jump = iters*iters;
            slong iters2;
      
            for (i = 0; i < iters; i++)
            {
               qfb_pow_ui(pow, pow, n, jump);
               qfb_pow_ui(twopow, pow, n, 4096);
               if (qfb_is_principal_form(twopow, n)) /* found factor */
                  break;
               iters2 = qfb_hash_find(qhash, twopow, depth); 
               if (iters2 != -1) /* found factor */
               {
                  if (fmpz_sgn(qhash[iters2].q->b) == fmpz_sgn(twopow->b))
                     qfb_inverse(qhash[iters2].q2, qhash[iters2].q2);

                  qfb_nucomp(pow, pow, qhash[iters2].q2, n, L);
                  qfb_reduce(pow, pow, n);

                  break;
               }
            }
            
            if (i == iters)
               done = 0;
         }

done:
         if (done)
         {
            for (i = 0; i < 12; i++)
            {
               qfb_pow_ui(twopow, pow, n, 2);
               if (qfb_is_principal_form(twopow, n))
               {
                  if (fmpz_cmpabs(pow->a, pow->b) != 0)
                  {
                     fmpz_abs(pow->b, pow->b);
                     fmpz_sub(g, pow->b, pow->a);
                     fmpz_sub(pow->a, g, pow->a);
                  }

                  fmpz_gcd(g, pow->a, n0);

                  if (!fmpz_is_one(g)) /* Success! */
                  {
                     printf("Factor: ");
                     fmpz_print(g);
                     printf("\n");
                     return 0;
                  }

                  done = 0;
                  break;
               }
               qfb_set(pow, twopow);
            }
         } 
   
         mult += 2;

         qfb_hash_clear(qhash, depth);

      }

      iters *= 2;
      jmax = (slong) (jmax * 1.15);

      mult = iters/10;
      mult |= 1;
   }

   qfb_clear(pow);
   qfb_clear(twopow);
   
   fmpz_clear(g);
   fmpz_clear(n);
   fmpz_clear(n0);
   fmpz_clear(p);
   fmpz_clear(L);

   return 0;
}
