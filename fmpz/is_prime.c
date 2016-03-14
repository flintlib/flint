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

    Copyright (C) 2014 William Hart
   
******************************************************************************/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "aprcl.h"

int fmpz_is_prime(const fmpz_t n)
{
   double logd = log(fmpz_get_d(n));
   ulong p, ppi, limit = (ulong) (logd*logd*logd/100.0) + 20;
   ulong * pp1, * pm1;
   slong i, l, num, num_pp1, num_pm1;
   const ulong * primes; 
   const double * pinv;

   fmpz_t F1, Fsqr, Fcub, R;
   int res = -1;
   
   if (fmpz_cmp_ui(n, 1) <= 0)
      return 0;
   
   if (fmpz_is_even(n))
      return (fmpz_cmp_ui(n, 2) == 0);

   if (fmpz_is_square(n))
      return 0;

   fmpz_init(F1);
   fmpz_init(R);
   fmpz_init(Fsqr);
   fmpz_init(Fcub);

   for (l = 0; l < 4 && res == -1; l++, limit *= 10)
   {
      num_pm1 = num_pp1 = 0;

      /* number of primes multiplied that will fit in a word */
      num = FLINT_BITS/FLINT_BIT_COUNT(limit);

      /* compute remainders of n mod p for primes p up to limit (approx.) */

      n_prime_pi_bounds(&ppi, &ppi, limit); /* precompute primes */
      primes = n_primes_arr_readonly(ppi + FLINT_BITS);
      pinv = n_prime_inverses_arr_readonly(ppi + FLINT_BITS);
   
      pm1 = _nmod_vec_init(2 + (ulong) logd); /* space for primes dividing n - 1 */
      pp1 = _nmod_vec_init(2 + (ulong) logd); /* space for primes dividing n + 1 */

      while (primes[0] < limit)
      {
         /* multiply batch of primes */
      
         p = primes[0];
         for (i = 1; i < num; i++)
            p *= primes[i];

         /* multi-modular reduction */

         p = fmpz_tdiv_ui(n, p);

         /* check for factors */
         for (i = 0; i < num; i++)
         {
            ulong r = n_mod2_precomp(p, primes[i], pinv[i]);

            if (r == 1) /* n - 1 = 0 mod p */
               pm1[num_pm1++] = primes[i];

            if (r == primes[i] - 1) /* n + 1 = 0 mod p */
               pp1[num_pp1++] = primes[i];
         }

         /* get next batch of primes */
         primes += num;
         pinv += num;
      }

      /* p - 1 test */
      res = fmpz_is_prime_pocklington(F1, R, n, pm1, num_pm1);

      if (res == 1)
      {
         fmpz_mul(Fsqr, F1, F1);
         if (fmpz_cmp(Fsqr, n) < 0)
         {
            fmpz_mul(Fcub, Fsqr, F1);
            if (fmpz_cmp(Fcub, n) >= 0) /* Brillhart, Lehmer, Selfridge test */
            {
               fmpz_t n1, c2, c1;

               fmpz_init(n1);
               fmpz_init(c2);
               fmpz_init(c1);

               fmpz_sub_ui(n1, n, 1); /* n is 1 mod F1 */
               fmpz_tdiv_q(n1, n1, F1);

               fmpz_tdiv_qr(c2, c1, n1, F1); /* Let n = c2*F^2 + c1*F + 1 */

               fmpz_mul(c1, c1, c1); /* check if c1^2 - 4*c2 is a square */
               fmpz_submul_ui(c1, c2, 4);

               if (fmpz_is_square(c1))
                  res = 0;
               /* else n is prime (res == 1) */

               fmpz_clear(n1);
               fmpz_clear(c2);
               fmpz_clear(c1);
            } else /* p + 1 test */
            {
               fmpz_t F2, Fm1;
         
               fmpz_init(F2);
               fmpz_init(Fm1);
               
               res = fmpz_is_prime_morrison(F2, R, n, pp1, num_pp1);

               if (res == 1)
               {
                  fmpz_sub_ui(Fm1, F2, 1); /* need F2 - 1 > sqrt(n) */
                  fmpz_mul(Fsqr, Fm1, Fm1);

                  if (fmpz_cmp(Fsqr, n) <= 0)
                  {
                     fmpz_mul(Fcub, Fsqr, Fm1);

                     if (fmpz_cmp(Fcub, n) > 0) /* Improved n + 1 test */
                     {
                        fmpz_t r1, r0, b, r, t;

                        fmpz_init(r1);
                        fmpz_init(r0);
                        fmpz_init(b);
                        fmpz_init(r);
                        fmpz_init(t);

                        fmpz_tdiv_qr(r1, r0, R, F2); /* R = r1*F2 + r0 */

                        /* check if x^2 + r0*x - r1 has positive integral root */
                        fmpz_mul(t, r0, r0); /* b = sqrt(r0^2 - 4(-r1)) */
                        fmpz_addmul_ui(t, r1, 4);
                        fmpz_sqrtrem(b, r, t);

                        if (fmpz_is_zero(r) && fmpz_cmp(b, r0) > 0) /* if so, composite */
                           res = 0;

                        /* check if x^2 + (r0 - F2)*x - r1 - 1 has positive integral root */
                        fmpz_sub(r0, r0, F2);
                        fmpz_add_ui(r1, r1, 1);
                        
                        fmpz_mul(t, r0, r0); /* b = sqrt((r0 - F2)^2 - 4(-r1 - 1)) */
                        fmpz_addmul_ui(t, r1, 4);
                        fmpz_sqrtrem(b, r, t);

                        if (fmpz_is_zero(r) && fmpz_cmp(b, r0) > 0) /* if so, composite */
                           res = 0;

                        fmpz_clear(t);
                        fmpz_clear(b);
                        fmpz_clear(r);
                        fmpz_clear(r1);
                        fmpz_clear(r0);
                     } else /* Brillhart, Lehmer, Selfridge combined p-1, p+1 test */
                     {
                        fmpz_t F, nmodF;
                        
                        fmpz_init(F);
                        
                        fmpz_mul(F, F1, F2); /* F = lcm(F1, F2), F1 | n - 1, F2 | n + 1 */
                        if (fmpz_is_even(F1) && fmpz_is_even(F2))
                           fmpz_tdiv_q_2exp(F, F, 1);

                        fmpz_mul(Fsqr, F, F);

                        if (fmpz_cmp(Fsqr, n) > 0) /* lcm(F1, F2) > sqrt(n) */
                        {
                            fmpz_init(nmodF);
                          
                            fmpz_mod(nmodF, n, F); /* check n mod F not factor of n */
                            
                            if (!fmpz_equal(nmodF, n) && !fmpz_is_one(nmodF) 
                              && fmpz_divisible(n, nmodF))
                               res = 0;

                            fmpz_clear(nmodF);
                        } else
                        {
                           fmpz_t d;
                           
                           fmpz_init(d);
                              
                           fmpz_mul(Fcub, Fsqr, F);
                           
                           if (fmpz_cmp(Fcub, n) > 0) /* Lenstra's divisors in residue class */
                           {
                              fmpz_t r;

                              fmpz_init(r);
                              
                              fmpz_set_ui(r, 1);
                              if (fmpz_divisor_in_residue_class_lenstra(d, n, r, F))
                                 res = 0;

                              fmpz_mod(r, n, F);
                              if (fmpz_divisor_in_residue_class_lenstra(d, n, r, F))
                                 res = 0;

                              fmpz_clear(r);
                           } else /* apr-cl primality test */
                           {
                              res = is_prime_aprcl(n);
                           }

                           fmpz_clear(d);
                        }

                        fmpz_clear(F);
                     }
                  }
                  /* else n is prime, i.e. res = 1 */
               }

               fmpz_clear(F2);
               fmpz_clear(Fm1);
            }
         }
      } 

      _nmod_vec_clear(pm1);
      _nmod_vec_clear(pp1);

   }

   fmpz_clear(F1);
   fmpz_clear(R);
   fmpz_clear(Fsqr);
   fmpz_clear(Fcub);

   return res;    
}
