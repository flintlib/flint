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

#define NUM_POLYS 7 /* x^3 - 1, x^4 - 1, x^6 - 1, x^5 - 1, x^8 - 1, x^12 - 1, x^10 - 1 */

static slong pol_degrees[NUM_POLYS] = { 3, 4, 6, 5, 8, 12, 10 };

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

void __print_fmpz_list(fmpz * arr, slong len)
{
   slong i;
   
   printf("[");
   for (i = 0; i < len - 1; i++)
      fmpz_print(arr + i), printf(", ");
   if (len)
      fmpz_print(arr + i);
   printf("]\n");
}

void __combine_residues(fmpz * r, fmpz * r1, fmpz_t F1, slong r1n, 
                                  fmpz * r2, fmpz_t F2, slong r2n)
{
   slong i, j;
   fmpz * t = _fmpz_vec_init(r1n*r2n);

   for (i = 0; i < r1n; i++)
   {
      for (j = 0; j < r2n; j++)
      {
         if (fmpz_is_zero(r2 + j))
            fmpz_set(t + i*r2n + j, r1 + i);
         else if (fmpz_is_zero(r1 + i))
            fmpz_set(t + i*r2n + j, r2 + j);
         else
            fmpz_CRT(t + i*r2n + j, r1 + i, F1, r2 + j, F2, 0);
      }
   }

   for (i = 0; i < r1n*r2n; i++)
      fmpz_set(r + i, t + i);

   _fmpz_vec_clear(t, r1n*r2n);
}

void __combine_residues2(fmpz * r, fmpz * r1, fmpz_t F1, slong r1n, 
                                   fmpz * r2, fmpz_t F2, slong r2n)
{
   slong i, j;

   for (i = 0, j = 0; i < r1n; i++, j++)
   {
      if (j == r2n) /* repeat same residues */
         j = 0;
      
      if (fmpz_is_zero(r2 + j))
         fmpz_set(r + i, r1 + i);
      else if (fmpz_is_zero(r1 + i))
         fmpz_set(r + i, r2 + j);
      else
         fmpz_CRT(r + i, r1 + i, F1, r2 + j, F2, 0);
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
   double logd = log(fmpz_get_d(n));
   ulong p, ppi, limit = (ulong) (logd*logd*logd/100.0) + 20;
   ulong * pp1, * pm1, ** pk1, * pk1n;
   slong i, j, l, num, num_pp1, num_pm1;
   const ulong * primes; 
   const double * pinv;

   fmpz_t F1, Fsqr, Fcub, R;
   int num_polys = 3;
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

   for (l = 0; l < 4 && res == -1; l++, limit *= 10, 
      num_polys = FLINT_MIN(num_polys + 2, NUM_POLYS))
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

      /* allocate space for prime arrays */
      pk1 = flint_malloc(num_polys*sizeof(ulong *));
      pk1n = flint_malloc(num_polys*sizeof(ulong));
      
      for (j = 0; j < num_polys; j++)
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
            ulong rk, r2, r3, r4, r = n_mod2_precomp(p, primes[i], pinv[i]);

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

            if (num_polys > 3)
            {
               /* check for n^4 + n^3 + n^2 + n + 1 = 0 mod p */
               r4 = n_mulmod_precomp(r2, r2, primes[i], pinv[i]);
               r3 = n_mulmod_precomp(r2, r, primes[i], pinv[i]);

               rk = n_addmod(r4, r3, primes[i]);
               rk = n_addmod(rk, r2, primes[i]);
               rk = n_addmod(rk, r, primes[i]);

               if (rk == primes[i] - 1) /* n^4 + n^3 + n^2 + n + 1 = 0 mod p */
                  pk1[3][pk1n[3]++] = primes[i];  

               /* check for n^4 + 1 = 0 mod p */
               if (r4 == primes[i] - 1) /* n^4 + 1 = 0 mod p */
                  pk1[4][pk1n[4]++] = primes[i];
            }

            if (num_polys > 5)
            {
               /* check for n^4 - n^2 + 1 = 0 mod p */
               rk = n_submod(r4, r2, primes[i]);
            
               if (rk == primes[i] - 1) /* n^4 - n^2 + 1 = 0 mod p */
                  pk1[5][pk1n[5]++] = primes[i];
            
               /* check for n^4 - n^3 + n^2 - n + 1 = 0 mod p */
               rk = n_submod(r4, r3, primes[i]);
               rk = n_addmod(rk, r2, primes[i]);
               rk = n_submod(rk, r, primes[i]);

               if (rk == primes[i] - 1) /* n^4 - n^3 + n^2 - n + 1 = 0 mod p */
                  pk1[6][pk1n[6]++] = primes[i];
            }
         }

         /* get next batch of primes */
         primes += num;
         pinv += num;
      }

      /* remove primes already dealt with at each size k */      
      for (j = 0; j < num_polys; j++)
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
                              fmpz * Fk;
                              fmpz ** rk = (fmpz **) flint_malloc((num_polys + 1)*sizeof(fmpz *));
                              fmpz ** racc = (fmpz **) flint_malloc((num_polys + 1)*sizeof(fmpz *));
                              slong * rkn = (slong *) flint_malloc((num_polys + 1)*sizeof(slong));
                              slong num_residues;
                              
                              Fk = _fmpz_vec_init(num_polys + 1);

                              rk[0] = _fmpz_vec_init(2);
                              for (i = 0; i < num_polys; i++)
                                 rk[i + 1] = _fmpz_vec_init(pol_degrees[i]);
                              
                              if (fmpz_is_one(F))
                                 fmpz_set_ui(rk[0] + 0, 0);
                              else
                                 fmpz_set_ui(rk[0] + 0, 1);
                              fmpz_mod(rk[0] + 1, n, F);
                              rkn[0] = fmpz_is_one(rk[0] + 1) || fmpz_is_one(F) ? 1 : 2;

                              fmpz_set(Fk + 0, F);

                              res = -1;

                              for (j = 0; res == -1 && j < num_polys; j++)
                              {
                                 if ((pol_degrees[j]) % 2 == 0)
                                    num_residues = pol_degrees[j];
                                 else
                                    num_residues = pol_degrees[j]*rkn[0];

                                 for (i = 0; i < j; i++)
                                 {
                                    if ((pol_degrees[j] % pol_degrees[i]) != 0)
                                       num_residues *= pol_degrees[i + 1];
                                 }
                                    
                                 racc[j + 1] = _fmpz_vec_init(num_residues);                                 
                                 rkn[j + 1] = num_residues;
                              
                                 res = 1;
                                 
                                 res = fmpz_is_prime_lenstra(Fk + j + 1, rk[j + 1], n, pk1[j], pk1n[j], pol_degrees[j]);

                                 if ((pol_degrees[j] % 2) == 0)
                                 {
                                    __combine_residues2(racc[j + 1], rk[j + 1], Fk + j + 1, pol_degrees[j], rk[0], Fk + 0, rkn[0]);
                                    fmpz_mul(F, Fk + j + 1, Fk + 0);
                                 } else
                                 {
                                    __set_residues(racc[j + 1], rk[j + 1], pol_degrees[j]);
                                    fmpz_set(F, Fk + j + 1);
                                 }

                                 for (i = 0; i < j; i++)
                                 {
                                    if ((pol_degrees[j] % pol_degrees[i]) == 0)
                                    {
                                       __combine_residues2(racc[j + 1], racc[j + 1], F, pol_degrees[j], rk[i + 1], Fk + i + 1, pol_degrees[i]);
                                       fmpz_mul(F, F, Fk + i + 1);
                                    }
                                 }

                                 num_residues = pol_degrees[j];
                                 
                                 if ((pol_degrees[j] % 2) != 0)
                                 {
                                    __combine_residues(racc[j + 1], racc[j + 1], F, num_residues, rk[0], Fk + 0, rkn[0]);
                                    fmpz_mul(F, F, Fk + 0);
                                    if (!fmpz_is_one(Fk + 0))
                                       num_residues *= rkn[0];
                                 } 

                                 for (i = 0; i < j; i++)
                                 {
                                    if ((pol_degrees[j] % pol_degrees[i]) != 0)
                                    {
                                       __combine_residues(racc[j + 1], racc[j + 1], F, num_residues, rk[i + 1], Fk + i + 1, pol_degrees[i]);
                                       fmpz_mul(F, F, Fk + i + 1);
                                       if (!fmpz_is_one(Fk + i + 1))
                                          num_residues *= pol_degrees[i];
                                    }
                                 }

                                 fmpz_mul(Fsqr, F, F);

                                 if (fmpz_cmp(Fsqr, n) >= 0) 
                                 {
                                    for (i = 0; i < num_residues; i++)
                                    {
                                       if (!fmpz_is_one(racc[j + 1] + i) && !fmpz_equal(racc[j + 1] + i, n) 
                                         && fmpz_divisible(n, racc[j + 1] + i))
                                         res = 0;
                                    }
                                 } else
                                 {
                                    fmpz_mul(Fcub, Fsqr, F);
                                    
                                    if (fmpz_cmp(Fcub, n) >= 0)
                                    {
                                       for (i = 0; i < num_residues; i++)
                                       {
                                          if (fmpz_divisor_in_residue_class_lenstra(d, n, racc[j + 1] + i, F))
                                             res = 0;
                                       }
                                    } else
                                       res = -1;
                                 }
                              }
                              
                              for (i = 0; i < j; i++)
                                 _fmpz_vec_clear(racc[i + 1], rkn[i + 1]);
                              
                              _fmpz_vec_clear(rk[0], 2);
                              for (i = 0; i < num_polys; i++)
                                 _fmpz_vec_clear(rk[i + 1], pol_degrees[i]);
                              
                              _fmpz_vec_clear(Fk, num_polys + 1);

                              flint_free(rk);
                              flint_free(racc);
                              flint_free(rkn);
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
      for (j = 0; j < num_polys; j++)
         _nmod_vec_clear(pk1[j]);

      flint_free(pk1);
      flint_free(pk1n);

      _nmod_vec_clear(pm1);
      _nmod_vec_clear(pp1);

   }

   fmpz_clear(F1);
   fmpz_clear(R);
   fmpz_clear(Fsqr);
   fmpz_clear(Fcub);

   return res;
      
}

#undef NUM_POLYS