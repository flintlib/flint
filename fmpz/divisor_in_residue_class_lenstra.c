/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

int fmpz_divisor_in_residue_class_lenstra(fmpz_t fac, const fmpz_t n, const fmpz_t r, const fmpz_t s)
{
   fmpz_t r1, r2, a0, a1, b0, b1, c0, c1, q, d, d1, s1, s2, ns2;
   slong i;
   int res = 0;

   fmpz_init(r1);
   fmpz_init(r2);
   fmpz_init(a0);
   fmpz_init(a1);
   fmpz_init(b0);
   fmpz_init(b1);
   fmpz_init(c0);
   fmpz_init(c1);
   fmpz_init(q);
   fmpz_init(d);
   fmpz_init(d1);
   fmpz_init(s1);
   fmpz_init(s2);
   fmpz_init(ns2);

   /* ns2 = n/s^2 */
   fmpz_mul(ns2, s, s);
   fmpz_tdiv_q(ns2, n, ns2);

   /* initialise */
   fmpz_invmod(r1, r, s); /* r1 = r^-1 mod s */
   fmpz_mul(r2, r1, n);
   fmpz_mod(r2, r2, s); /* r2 = r1*n mod s */
   
   fmpz_set(a0, s); /* a0 = s */
   fmpz_mul(a1, r1, r2); /* a1 = r1*r2 mod s */
   fmpz_mod(a1, a1, s);
   if (fmpz_is_zero(a1))
      fmpz_add(a1, a1, s); /* 0 < a1 <= s */
   
   fmpz_zero(b0);
   fmpz_one(b1);
   
   fmpz_zero(c0);
   fmpz_mul(c1, r, r2); /* c1 = (n - r*r2)/s * r1 mod s */
   fmpz_sub(c1, n, c1);
   fmpz_divexact(c1, c1, s);
   fmpz_mul(c1, c1, r1);
   fmpz_mod(c1, c1, s);


   /* deal with a0, b0, c0 */
   if (!fmpz_is_one(r) && !fmpz_equal(n, r) && fmpz_divisible(n, r))
   {
      fmpz_set(fac, r);
      res = 1;
      goto cleanup;
   }

   for (i = 1; ; i++)
   {
      if ((i & 1) == 0)
      {
         fmpz_mod(s1, c1, s);
         fmpz_neg(s2, s); 
      } else
      {
         fmpz_mul(s2, a1, b1); /* s2 = a1*b1 */
         fmpz_add(s1, s2, ns2); /* s1 = a1*b1 + n/s^2 */

         fmpz_mod(q, s1, s); /* s1 largest integer < a1*b1 + n/s^2 congruent to c1 mod s */
         if (fmpz_cmp(q, c1) < 0)
            fmpz_sub(s1, s1, s);
         fmpz_sub(s1, s1, q);
         fmpz_add(s1, s1, c1);

         fmpz_add(s2, s2, s2); /* s2 = 2*a1*b1 - 1 */
         fmpz_sub_ui(s2, s2, 1);
      }

      while (fmpz_cmp(s1, s2) > 0) /* for each value s1 in range */
      {
         fmpz_mul(d, s, s1); /* d = (s1*s + a1*r + b1*r2)^2 - 4*a1*b1*n */
         fmpz_addmul(d, a1, r);
         fmpz_addmul(d, b1, r2);
         fmpz_set(d1, d); /* d1 = s1*s + a1*r + b1*r2 */
         fmpz_mul(d, d, d);
         fmpz_mul(q, a1, b1);
         fmpz_mul(q, q, n);
         fmpz_submul_ui(d, q, 4);

         if (fmpz_is_square(d)) /* divisor exists, roots are (d1 +/- sqrt(d))/2 */
         {
            fmpz_sqrt(d, d);
            fmpz_add(d, d, d1);
            fmpz_tdiv_q_2exp(d, d, 1);

            if (fmpz_is_zero(a1))
            {
               fmpz_tdiv_q(fac, s1, b1); /* y*b1 = s1 */
               fmpz_mul(fac, fac, s); /* check if ys + r2 is factor */
               fmpz_add(fac, fac, r2);
               fmpz_abs(fac, fac);

               if (!fmpz_is_zero(fac) && !fmpz_is_one(fac) && !fmpz_equal(fac, n) && fmpz_divisible(n, fac))
               {
                  res = 1;
                  break;
               }
            } else if (fmpz_is_zero(b1))
            {
               fmpz_tdiv_q(fac, s1, a1); /* x*a1 = s1 */
               fmpz_mul(fac, fac, s); /* check if xs + r is factor */
               fmpz_add(fac, fac, r);
               fmpz_abs(fac, fac);

               if (!fmpz_is_zero(fac) && !fmpz_is_one(fac) && !fmpz_equal(fac, n) && fmpz_divisible(n, fac))
               {
                  res = 1;
                  break;
               }
            } else
            {
               /* either d/a1 or d/b1 is a divisor of n */
               fmpz_tdiv_q(fac, d, a1);
               fmpz_abs(fac, fac);
               if (!fmpz_is_zero(fac) && !fmpz_is_one(fac) && !fmpz_equal(fac, n) && fmpz_divisible(n, fac))
               {
                  res = 1;
                  break;
               }

               fmpz_tdiv_q(fac, d, b1);
               fmpz_abs(fac, fac);
               if (!fmpz_is_zero(fac) && !fmpz_is_one(fac) && !fmpz_equal(fac, n) && fmpz_divisible(n, fac))
               {
                  res = 1;
                  break;
               }
            }
         }

         fmpz_sub(s1, s1, s);
      }
      
      if (fmpz_is_zero(a1) || res == 1) /* Euclidean chain has terminated */
         break;
      
      /* 
         Euclidean chain:
         
         a1, a0 = a0 - q*a1, a1, 
         where 0 <= a1 < a0 for i even
         and 0 < a1 <= a0 for i odd
      */
      fmpz_tdiv_qr(q, a0, a0, a1);
      if ((i & 1) == 0 && fmpz_is_zero(a0))
      {
         fmpz_sub_ui(q, q, 1);
         fmpz_add(a0, a0, a1);
      }
      fmpz_swap(a0, a1);

      /* b1, b0 = b0 - q*b1, b1 */
      fmpz_submul(b0, q, b1);
      fmpz_swap(b0, b1);

      /* c1, c0 = c0 - q*c1, c1 mod s */
      fmpz_submul(c0, q, c1);
      fmpz_mod(c0, c0, s);
      fmpz_swap(c0, c1);
   }

cleanup:

   fmpz_clear(r1);
   fmpz_clear(r2);
   fmpz_clear(a0);
   fmpz_clear(a1);
   fmpz_clear(b0);
   fmpz_clear(b1);
   fmpz_clear(c0);
   fmpz_clear(c1);
   fmpz_clear(q);
   fmpz_clear(d);
   fmpz_clear(d1);
   fmpz_clear(s1);
   fmpz_clear(s2);
   fmpz_clear(ns2);

   return res;
}
