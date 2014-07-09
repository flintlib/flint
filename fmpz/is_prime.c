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

#define NUM_POLYS 4 /* x^3 - 1, x^4 - 1, x^6 - 1, x^5 - 1 */

static slong pol_degrees[4] = { 3, 4, 6, 5 };

slong __primes_remove(ulong * arr1, slong len1, ulong * arr2, slong len2)
{
   slong i, j, k;

   for (i = 0, j = 0, k = 0; i < len1; )
   {
      if (k >= len2 || arr1[i] < arr2[k])
         arr1[j++] = arr1[i++];
      else if (arr1[i] == arr2[k])
         i++, k++;
      else
         k++;
   }

   return j;
}

void __print_list(ulong * arr, slong len)
{
   slong i;
   
   printf("[");
   for (i = 0; i < len - 1; i++)
      printf("%ld, ", arr[i]);
   if (len)
      printf("%ld", arr[i]);
   printf("]\n");
}

void __combine_residues(fmpz * r, fmpz * r1, fmpz_t F1, slong r1n, 
                                  fmpz * r2, fmpz_t F2, slong r2n)
{
   slong i, j;

   for (i = 0; i < r1n; i++)
   {
      for (j = 0; j < r2n; j++)
         fmpz_CRT(r + i*r2n + j, r1 + i, F1, r2 + j, F2, 0);
   }
}

void __set_residues(fmpz * r, fmpz * r1, slong r1n)
{
   slong i;

   for (i = 0; i < r1n; i++)
      fmpz_set(r + i, r1 + i);
}

int fmpz_is_prime(const fmpz_t n)
{
   int res = 0;
   
   if (fmpz_cmp_ui(n, 1) <= 0)
      return 0;
   
   if (fmpz_is_even(n))
      return (fmpz_cmp_ui(n, 2) == 0);

   if (fmpz_is_square(n))
      return 0;

   {
      double logd = log(fmpz_get_d(n));
      ulong p, ppi, limit = (ulong) (logd*logd*logd/10.0) + 2;
      ulong * pp1, * pm1, ** pk1, * pk1n;
      slong i, j, num, num_pp1 = 0, num_pm1 = 0;
      const ulong * primes; 
      const double * pinv;

      fmpz_t F1, Fsqr, Fcub, R;

      fmpz_init(F1);
      fmpz_init(R);
      fmpz_init(Fsqr);
      fmpz_init(Fcub);

      /* number of primes multiplied that will fit in a word */
      num = FLINT_BITS/FLINT_BIT_COUNT(limit);

      /* compute remainders of n mod p for primes p up to limit (approx.) */

      n_prime_pi_bounds(&ppi, &ppi, limit); /* precompute primes */
      primes = n_primes_arr_readonly(ppi + FLINT_BITS);
      pinv = n_prime_inverses_arr_readonly(ppi + FLINT_BITS);
   
      pm1 = _nmod_vec_init(2 + (ulong) logd); /* space for primes dividing n - 1 */
      pp1 = _nmod_vec_init(2 + (ulong) logd); /* space for primes dividing n + 1 */

      /* allocate space for prime arrays */
      pk1 = flint_malloc(NUM_POLYS*sizeof(ulong *));
      pk1n = flint_malloc(NUM_POLYS*sizeof(ulong));
      
      for (j = 0; j < NUM_POLYS; j++)
      {
         pk1[j] = _nmod_vec_init(pol_degrees[j] * (ulong) logd);
         pk1n[j] = 0;
      }

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
            ulong rk, r2, r4, r = n_mod2_precomp(p, primes[i], pinv[i]);

            if (r == 1) /* n - 1 = 0 mod p */
               pm1[num_pm1++] = primes[i];

            if (r == primes[i] - 1) /* n + 1 = 0 mod p */
               pp1[num_pp1++] = primes[i];

            /* n^2 mod p */
            r2 = n_mulmod_precomp(r, r, primes[i], pinv[i]);
               
            /* check for n^2 + n + 1 = 0 mod p */
            rk = n_addmod(r2, r, primes[i]);

            if (rk == primes[i] - 1) /* n^2 + n + 1 = 0 mod p */
               pk1[0][pk1n[0]++] = primes[i];

            /* check for n^2 + 1 = 0 mod p */
            if (r2 == primes[i] - 1) /* n^2 + 1 = 0 mod p */
               pk1[1][pk1n[1]++] = primes[i];

            /* check for n^2 - n + 1 = 0 mod p */
            rk = n_submod(r2, r, primes[i]);
            
            if (rk == primes[i] - 1) /* n^2 - n + 1 = 0 mod p */
               pk1[2][pk1n[2]++] = primes[i];  

            /* check for n^4 + 1 = 0 mod p */
            r4 = n_mulmod_precomp(r2, r2, primes[i], pinv[i]);

            if (r4 == primes[i] - 1) /* n^4 + 1 = 0 mod p */
               pk1[3][pk1n[3]++] = primes[i];  
         }

         /* get next batch of primes */
         primes += num;
         pinv += num;
      }

      /* remove primes already dealt with at each size k */      
      for (j = 0; j < NUM_POLYS; j++)
      {
          pk1n[j] = __primes_remove(pk1[j], pk1n[j], pm1, num_pm1);
          pk1n[j] = __primes_remove(pk1[j], pk1n[j], pp1, num_pp1);

          for (i = 0; i < j; i++)
             pk1n[j] = __primes_remove(pk1[j], pk1n[j], pk1[i], pk1n[i]);
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
                           } else /* Lenstra finite fields primality test */
                           {
                              fmpz_t Fk;
                              fmpz * r2, * r3, * r4, * r6, * r5, * r23, * r234, * r2346, * r23465;
                              slong r2n, r3n, r4n, r6n, r5n, r23n, r234n, r2346n, r23465n;

                              fmpz_init(Fk);
                              
                              r2 = _fmpz_vec_init(2);
                              r3 = _fmpz_vec_init(3);
                              r4 = _fmpz_vec_init(4);
                              r6 = _fmpz_vec_init(6);
                              r5 = _fmpz_vec_init(5);
                              
                              fmpz_set_ui(r2 + 0, 1);
                              fmpz_mod(r2 + 1, n, F);
                              r2n = fmpz_is_one(r2 + 1) ? 1 : 2;

                              res = -1;

                              /* p^2 + p + 1 test */
                              r23 = _fmpz_vec_init(3*r2n);                                 
                              
                              if (pk1n[0] != 0)
                              {
                                 res = 1;
                                 
                                 r23n = 3*r2n;

                                 res = fmpz_is_prime_lenstra3(Fk, r3, n, pk1[0], pk1n[0]);
                                 r3n = 3;

                                 if (!fmpz_is_one(Fk))
                                 {
                                    __combine_residues(r23, r3, Fk, r3n, r2, F, r2n);

                                    fmpz_mul(F, F, Fk);
                                    fmpz_mul(Fsqr, F, F);

                                    if (fmpz_cmp(Fsqr, n) >= 0) 
                                    {
                                       for (i = 0; i < r23n; i++)
                                       {
                                          if (!fmpz_is_one(r23 + i) && !fmpz_equal(r23 + i, n) 
                                            && fmpz_divisible(n, r23 + i))
                                            res = 0;
                                       }
                                    } else
                                    {
                                       fmpz_mul(Fcub, Fsqr, F);
                                    
                                       if (fmpz_cmp(Fcub, n) >= 0)
                                       {
                                          for (i = 0; i < r23n; i++)
                                          {
                                             if (fmpz_divisor_in_residue_class_lenstra(d, n, r23 + i, F))
                                                res = 0;
                                          }
                                       } else
                                          res = -1;
                                    }
                                 } else
                                 {
                                    __set_residues(r23, r2, r2n);
                                    r23n = r2n;
                                 }
                              } else
                              {
                                 __set_residues(r23, r2, r2n);
                                 r23n = r2n;
                              }

                              /* p^2 + 1 test */
                              r234 = _fmpz_vec_init(4*r23n);

                              if (res == -1 && pk1n[1] != 0)
                              {
                                 res = 1;
                                 
                                 r234n = 4*r23n;

                                 res = fmpz_is_prime_lenstra4(Fk, r4, n, pk1[1], pk1n[1]);
                                 r4n = 4;

                                 if (!fmpz_is_one(Fk))
                                 {
                                    __combine_residues(r234, r4, Fk, r4n, r23, F, r23n);

                                    fmpz_mul(F, F, Fk);
                                    fmpz_mul(Fsqr, F, F);

                                    if (fmpz_cmp(Fsqr, n) >= 0) 
                                    {
                                       for (i = 0; i < r234n; i++)
                                       {
                                          if (!fmpz_is_one(r234 + i) && !fmpz_equal(r234 + i, n) 
                                            && fmpz_divisible(n, r234 + i))
                                            res = 0;
                                       }
                                    } else
                                    {
                                       fmpz_mul(Fcub, Fsqr, F);
                                    
                                       if (fmpz_cmp(Fcub, n) >= 0)
                                       {
                                          for (i = 0; i < r234n; i++)
                                          {
                                             if (fmpz_divisor_in_residue_class_lenstra(d, n, r234 + i, F))
                                                res = 0;
                                          }
                                       } else
                                          res = -1;
                                    }
                                 } else
                                 {
                                    __set_residues(r234, r23, r23n);
                                    r234n = r23n;
                                 }
                              } else
                              {
                                 __set_residues(r234, r23, r23n);
                                 r234n = r23n;
                              }

                              /* p^2 - p + 1 test */
                              r2346 = _fmpz_vec_init(6*r234n);

                              if (res == -1 && pk1n[2] != 0)
                              {
                                 res = 1;
                                 
                                 r2346n = 6*r234n;

                                 res = fmpz_is_prime_lenstra6(Fk, r6, n, pk1[2], pk1n[2]);
                                 r6n = 6;

                                 if (!fmpz_is_one(Fk))
                                 {
                                    __combine_residues(r2346, r6, Fk, r6n, r234, F, r234n);

                                    fmpz_mul(F, F, Fk);
                                    fmpz_mul(Fsqr, F, F);

                                    if (fmpz_cmp(Fsqr, n) >= 0) 
                                    {
                                       for (i = 0; i < r2346n; i++)
                                       {
                                          if (!fmpz_is_one(r2346 + i) && !fmpz_equal(r2346 + i, n) 
                                            && fmpz_divisible(n, r2346 + i))
                                            res = 0;
                                       }
                                    } else
                                    {
                                       fmpz_mul(Fcub, Fsqr, F);
                                    
                                       if (fmpz_cmp(Fcub, n) >= 0)
                                       {
                                          for (i = 0; i < r2346n; i++)
                                          {
                                             if (fmpz_divisor_in_residue_class_lenstra(d, n, r2346 + i, F))
                                                res = 0;
                                          }
                                       } else
                                          res = -1;
                                    }
                                 } else
                                 {
                                    __set_residues(r2346, r234, r234n);
                                    r2346n = r234n;
                                 }
                              } else
                              {
                                 __set_residues(r2346, r234, r234n);
                                 r2346n = r234n;
                              }

                              /* p^4 + 1 test */
                              r23465 = _fmpz_vec_init(5*r2346n);

                              if (res == -1 && pk1n[2] != 0)
                              {
                                 res = 1;
                                 
                                 r23465n = 5*r2346n;

                                 res = fmpz_is_prime_lenstra5(Fk, r5, n, pk1[3], pk1n[3]);
                                 r5n = 5;

                                 if (!fmpz_is_one(Fk))
                                 {
                                    __combine_residues(r23465, r5, Fk, r5n, r2346, F, r2346n);

                                    fmpz_mul(F, F, Fk);
                                    fmpz_mul(Fsqr, F, F);

                                    if (fmpz_cmp(Fsqr, n) >= 0) 
                                    {
                                       for (i = 0; i < r23465n; i++)
                                       {
                                          if (!fmpz_is_one(r23465 + i) && !fmpz_equal(r23465 + i, n) 
                                            && fmpz_divisible(n, r23465 + i))
                                            res = 0;
                                       }
                                    } else
                                    {
                                       fmpz_mul(Fcub, Fsqr, F);
                                    
                                       if (fmpz_cmp(Fcub, n) >= 0)
                                       {
                                          for (i = 0; i < r23465n; i++)
                                          {
                                             if (fmpz_divisor_in_residue_class_lenstra(d, n, r23465 + i, F))
                                                res = 0;
                                          }
                                       } else
                                          res = -1;
                                    }
                                 } else
                                 {
                                    __set_residues(r23465, r2346, r2346n);
                                    r23465n = r2346n;
                                 }
                              } else
                              {
                                 __set_residues(r23465, r2346, r2346n);
                                 r23465n = r2346n;
                              }

                              _fmpz_vec_clear(r23, 3*r2n);
                              _fmpz_vec_clear(r234, 4*r23n);
                              _fmpz_vec_clear(r2346, 6*r234n);
                              _fmpz_vec_clear(r23465, 5*r2346n);

                              _fmpz_vec_clear(r3, 3);
                              _fmpz_vec_clear(r4, 4);
                              _fmpz_vec_clear(r6, 6);
                              _fmpz_vec_clear(r5, 5);

                              fmpz_clear(Fk);
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

      /* deallocate prime arrays */
      for (j = 0; j < NUM_POLYS; j++)
         _nmod_vec_clear(pk1[j]);

      flint_free(pk1);
      flint_free(pk1n);

      _nmod_vec_clear(pm1);
      _nmod_vec_clear(pp1);

      fmpz_clear(F1);
      fmpz_clear(R);
      fmpz_clear(Fsqr);
      fmpz_clear(Fcub);
   }

   return res;
      
}

#undef NUM_POLYS