/*
    Copyright (C) 2009, 2012, 2014 William Hart
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "mpn_extras.h"
#include "nmod_vec.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "aprcl.h"

int fmpz_is_prime(const fmpz_t n)
{
   double logd;
   ulong p, ppi, limit;
   ulong * pp1, * pm1;
   slong i, l, num, num_pp1, num_pm1;
   const ulong * primes;
   const double * pinv;

   fmpz_t F1, Fsqr, Fcub, R, t;
   int res = -1;

   if (fmpz_cmp_ui(n, 1) <= 0)
      return 0;

   if (fmpz_abs_fits_ui(n))
      return n_is_prime(fmpz_get_ui(n));

   if (fmpz_is_even(n))
      return 0;

   if (flint_mpn_factor_trial(COEFF_TO_PTR(*n)->_mp_d, COEFF_TO_PTR(*n)->_mp_size, 1, fmpz_bits(n)))
        return 0;

   /* todo: use fmpz_is_perfect_power? */
   if (fmpz_is_square(n))
      return 0;

   /* Fast deterministic Miller-Rabin test up to about 81 bits. This choice of
      bases certifies primality for n < 3317044064679887385961981;
      see https://doi.org/10.1090/mcom/3134 */
   fmpz_init(t);
   fmpz_tdiv_q_2exp(t, n, 64);
   if (fmpz_cmp_ui(t, 179817) < 0)
   {
      static const char bases[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 0 };

      for (i = 0; bases[i] != 0; i++)
      {
         fmpz_set_ui(t, bases[i]);
         if (!fmpz_is_strong_probabprime(n, t))
             return 0;  /* no need to clear t since it is small */
      }

      return 1;
   }

   /* Do a single base-2 test to rule out most composites */
   fmpz_set_ui(t, 2);
   if (!fmpz_is_strong_probabprime(n, t))
     return 0;

   fmpz_clear(t);

   logd = fmpz_dlog(n);
   limit = (ulong) (logd*logd*logd/100.0) + 20;

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
                              res = aprcl_is_prime(n);
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

   /* aprcl_is_prime() actually throws, but it does not hurt to have
      this fallback here */
   if (res < 0)
   {
      flint_printf("Exception in fmpz_is_prime: failed to prove ");
      fmpz_print(n);
      flint_printf(" prime or composite\n");
      flint_abort();
   }

   fmpz_clear(F1);
   fmpz_clear(R);
   fmpz_clear(Fsqr);
   fmpz_clear(Fcub);

   return res;
}

void _fmpz_np1_trial_factors(const fmpz_t n, mp_ptr pp1, slong * num_pp1, ulong limit)
{
   slong i, num;
   ulong ppi, p;
   const ulong * primes;
   const double * pinv;

   *num_pp1 = 0;

   /* number of primes multiplied that will fit in a word */

   num = FLINT_BITS/FLINT_BIT_COUNT(limit);

   /* compute remainders of n mod p for primes p up to limit (approx.) */

   n_prime_pi_bounds(&ppi, &ppi, limit); /* precompute primes */
   primes = n_primes_arr_readonly(ppi + FLINT_BITS);
   pinv = n_prime_inverses_arr_readonly(ppi + FLINT_BITS);

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

         if (r == primes[i] - 1) /* n + 1 = 0 mod p */
            pp1[(*num_pp1)++] = primes[i];
      }

      /* get next batch of primes */
      primes += num;
      pinv += num;
   }
}

int fmpz_is_prime_morrison(fmpz_t F, fmpz_t R, const fmpz_t n, mp_ptr pp1, slong num_pp1)
{
   slong i, d, bits;
   mp_limb_t a, b;
   fmpz_t g, q, r, ex, c, D, Dinv, A, B, Ukm, Ukm1, Um, Um1, Vm, Vm1, p;
   fmpz_factor_t fac;
   int res = 0, fac_found;

   fmpz_init(D);
   fmpz_init(Dinv);
   fmpz_init(A);
   fmpz_init(B);
   fmpz_init(p);
   fmpz_init(q);
   fmpz_init(r);
   fmpz_init(g);
   fmpz_init(c);
   fmpz_init(ex);
   fmpz_init(Um);
   fmpz_init(Um1);
   fmpz_init(Ukm);
   fmpz_init(Ukm1);
   fmpz_init(Vm);
   fmpz_init(Vm1);
   fmpz_factor_init(fac);

   fmpz_add_ui(R, n, 1); /* start with n + 1 */

   bits = fmpz_bits(R);

   for (i = 0; i < num_pp1; i++)
   {
      fmpz_set_ui(p, pp1[i]);
      d = fmpz_remove(R, R, p);
      _fmpz_factor_append_ui(fac, pp1[i], d);
   }

   srand(time(NULL));

   if (!fmpz_is_probabprime_BPSW(R))
   {
      if (bits > 150 && (fac_found = fmpz_factor_pp1(p, R, bits + 1000, bits/20 + 1000, rand()%100 + 3)
                    && fmpz_is_prime(p)))
      {
         d = fmpz_remove(R, R, p);
         _fmpz_factor_append(fac, p, d);

         if (fmpz_is_probabprime_BPSW(R)) /* fast test first */
         {
            if (fmpz_is_prime(R) == 1)
            {
               _fmpz_factor_append(fac, R, 1);
               fmpz_set_ui(R, 1);
            }
         }
      }
   } else
   {
      if (fmpz_is_prime(R) == 1)
      {
         _fmpz_factor_append(fac, R, 1);
         fmpz_set_ui(R, 1);
      }
   }

   /* compute product F of found primes */
   fmpz_set_ui(F, 1);
   for (i = 0; i < fac->num; i++)
   {
      if (fac->exp[i] == 1)
         fmpz_mul(F, F, fac->p + i);
      else
      {
         fmpz_pow_ui(ex, fac->p + i, fac->exp[i]);
         fmpz_mul(F, F, ex);
      }
   }

   /* Want D = A^2 - 4B where A = a, B = b, such that (D/n) = -1 */
   for (b = 1; ; b++)
   {
      fmpz_set_ui(B, b);
      fmpz_gcd(g, B, n);

      if (fmpz_equal(g, n)) /* need gcd(n, b) = 1 */
         continue;

      if (!fmpz_is_one(g)) /* found a factor of n */
      {
         res = 0;
         goto cleanup;
      }

      a = 2;
      do {
         a++;
         fmpz_set_ui(A, a);
         fmpz_mul_ui(D, A, a);
         fmpz_sub_ui(D, D, 4*b);
      } while (fmpz_jacobi(D, n) != -1);

      fmpz_invmod(Dinv, D, n);

      /* compute U((n+1)/F) mod n */
      fmpz_lucas_chain_full(Vm, Vm1, A, B, R, n);
      fmpz_lucas_chain_VtoU(Um, Um1, Vm, Vm1, A, B, Dinv, n);

      /* check U(n+1) = 0 mod n */
      fmpz_lucas_chain_mul(Ukm, Ukm1, Um, Um1, A, B, F, n);
      if (!fmpz_is_zero(Ukm))
      {
         res = 0;
         goto cleanup;
      }

      fmpz_set_ui(c, 1);

      /* find values U((n+1)/q) for each prime q dividing F */
      for (i = 0; i < fac->num; i++)
      {
         fmpz_tdiv_q(ex, F, fac->p + i);
         fmpz_lucas_chain_mul(Ukm, Ukm1, Um, Um1, A, B, ex, n);

         if (!fmpz_is_zero(Ukm))
         {
            fmpz_mul(c, c, Ukm);
            fmpz_mod(c, c, n);
         } else
            break;
      }

      if (i == fac->num) /* found valid base a */
         break;
   }

   /* check for factors of n */
   fmpz_gcd(g, n, c);
   res = fmpz_is_one(g);

cleanup:

   fmpz_factor_clear(fac);
   fmpz_clear(D);
   fmpz_clear(Dinv);
   fmpz_clear(A);
   fmpz_clear(B);
   fmpz_clear(c);
   fmpz_clear(ex);
   fmpz_clear(p);
   fmpz_clear(q);
   fmpz_clear(r);
   fmpz_clear(g);
   fmpz_clear(Um);
   fmpz_clear(Um1);
   fmpz_clear(Ukm);
   fmpz_clear(Ukm1);
   fmpz_clear(Vm);
   fmpz_clear(Vm1);

   return res;
}

void _fmpz_nm1_trial_factors(const fmpz_t n, mp_ptr pm1, slong * num_pm1, ulong limit)
{
   slong i, num;
   ulong ppi, p;
   const ulong * primes;
   const double * pinv;

   *num_pm1 = 0;

   /* number of primes multiplied that will fit in a word */

   num = FLINT_BITS/FLINT_BIT_COUNT(limit);

   /* compute remainders of n mod p for primes p up to limit (approx.) */

   n_prime_pi_bounds(&ppi, &ppi, limit); /* precompute primes */
   primes = n_primes_arr_readonly(ppi + FLINT_BITS);
   pinv = n_prime_inverses_arr_readonly(ppi + FLINT_BITS);

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
            pm1[(*num_pm1)++] = primes[i];
      }

      /* get next batch of primes */
      primes += num;
      pinv += num;
   }
}

int fmpz_is_prime_pocklington(fmpz_t F, fmpz_t R, const fmpz_t n, mp_ptr pm1, slong num_pm1)
{
   slong i, d, bits;
   ulong a;
   fmpz_t g, q, r, pow, pow2, ex, c, p;
   fmpz_factor_t fac;
   int res = 0, fac_found;

   fmpz_init(p);
   fmpz_init(q);
   fmpz_init(r);
   fmpz_init(g);
   fmpz_init(pow);
   fmpz_init(pow2);
   fmpz_init(c);
   fmpz_init(ex);
   fmpz_factor_init(fac);

   fmpz_sub_ui(R, n, 1); /* start with n - 1 */

   bits = fmpz_bits(R);

   for (i = 0; i < num_pm1; i++)
   {
      fmpz_set_ui(p, pm1[i]);
      d = fmpz_remove(R, R, p);
      _fmpz_factor_append_ui(fac, pm1[i], d);
   }

   srand(time(NULL));

   if (!fmpz_is_probabprime_BPSW(R))
   {
      if (bits > 150 && (fac_found = fmpz_factor_pp1(p, R, bits + 1000, bits/20 + 1000, rand()%100 + 3)
                    && fmpz_is_prime(p)))
      {
         d = fmpz_remove(R, R, p);
         _fmpz_factor_append(fac, p, d);

         if (fmpz_is_probabprime_BPSW(R)) /* fast test first */
         {
            if (fmpz_is_prime(R) == 1)
            {
               _fmpz_factor_append(fac, R, 1);
               fmpz_set_ui(R, 1);
            }
         }
      }
   } else
   {
      if (fmpz_is_prime(R) == 1)
      {
         _fmpz_factor_append(fac, R, 1);
         fmpz_set_ui(R, 1);
      }
   }

   /* compute product F of found primes */
   fmpz_set_ui(F, 1);
   for (i = 0; i < fac->num; i++)
   {
      if (fac->exp[i] == 1)
         fmpz_mul(F, F, fac->p + i);
      else
      {
         fmpz_pow_ui(pow, fac->p + i, fac->exp[i]);
         fmpz_mul(F, F, pow);
      }
   }

   for (a = 2; ; a++)
   {
      /* compute a^((n-1)/F) mod n */
      fmpz_set_ui(pow, a);
      fmpz_powm(pow, pow, R, n);

      /* check a^(n-1) = 1 mod n */
      fmpz_powm(pow2, pow, F, n);
      if (!fmpz_is_one(pow2))
      {
         res = 0;
         goto cleanup;
      }

      fmpz_set_ui(c, 1);

      /* find values a^((n-1)/q) - 1 for each prime q dividing F */
      for (i = 0; i < fac->num; i++)
      {
         fmpz_tdiv_q(ex, F, fac->p + i);
         fmpz_powm(pow2, pow, ex, n);
         fmpz_sub_ui(pow2, pow2, 1);
         if (fmpz_sgn(pow2) < 0)
            fmpz_add(pow2, pow2, n);

         if (!fmpz_is_zero(pow2))
         {
            fmpz_mul(c, c, pow2);
            fmpz_mod(c, c, n);
         } else
            break;
      }

      if (i == fac->num) /* found valid base a */
         break;
   }

   /* check for factors of n */
   fmpz_gcd(g, n, c);
   res = fmpz_is_one(g);

cleanup:

   fmpz_factor_clear(fac);
   fmpz_clear(pow);
   fmpz_clear(pow2);
   fmpz_clear(c);
   fmpz_clear(ex);
   fmpz_clear(p);
   fmpz_clear(q);
   fmpz_clear(r);
   fmpz_clear(g);

   return res;
}

#ifndef FLINT64
mp_limb_t flint_fmpz_pseudosquares[][3] =
{
   { 17, 0, 0 },
   { 73, 0, 0 },
   { 241, 0, 0 },
   { 1009, 0, 0 },
   { 2641, 0, 0 },
   { 8089, 0, 0 },
   { 18001, 0, 0 },
   { 53881, 0, 0 },
   { 87481, 0, 0 },
   { 117049, 0, 0 },
   { 515761, 0, 0 },
   { 1083289, 0, 0 },
   { 3206641, 0, 0 },
   { 3818929, 0, 0 },
   { 9257329, 0, 0 },
   { 22000801, 0, 0 },
   { 48473881, 0, 0 },
   { 48473881, 0, 0 },
   { 175244281, 0, 0 },
   { 427733329, 0, 0 },
   { 427733329, 0, 0 },
   { 898716289u, 0, 0 },
   { 2805544681u, 0, 0 },
   { 2805544681u, 0, 0 },
   { 2805544681u, 0, 0 },
   { 1720328849u, 2u, 0 },
   { 2141495009u, 5u, 0 },
   { 3553231785u, 19u, 0 },
   { 3553231785u, 19u, 0 },
   { 2991566689u, 45u, 0 },
   { 2991566689u, 45u, 0 },
   { 2804689073u, 668u, 0 },
   { 2804689073u, 668u, 0 },
   { 2804689073u, 668u, 0 },
   { 46910577u, 6112u, 0 },
   { 46910577u, 6112u, 0 },
   { 1079027281u, 26178u, 0 },
   { 1079027281u, 26178u, 0 },
   { 1079027281u, 26178u, 0 },
   { 3590018425u, 41661u, 0 },
   { 3590018425u, 41661u, 0 },
   { 2746102297u, 162087u, 0 },
   { 2746102297u, 162087u, 0 },
   { 1936779721u, 664710u, 0 },
   { 1070451441u, 1501768u, 0 },
   { 1070451441u, 1501768u, 0 },
   { 2061289617u, 2710474u, 0 },
   { 2061289617u, 2710474u, 0 },
   { 4235760785u, 44382509u, 0 },
   { 2312776601u, 45783875u, 0 },
   { 2678348049u, 165920782u, 0 },
   { 3315991761u, 413007985u, 0 },
   { 1567179849u, 541956877u, 0 },
   { 2273104657u, 1486621767u, 0 },
   { 3796117489u, 1867116582u, 0 },
   { 2425538377u, 2374430322u, 0 },
   { 2425538377u, 2374430322u, 0 },
   { 2425538377u, 2374430322u, 0 },
   { 3763487577u, 3377920039u, 3u },
   { 2972093785u, 1402148275u, 11u },
   { 2785759393u, 3968325740u, 28u },
   { 551239041u, 3335735663u, 50u },
   { 551239041u, 3335735663u, 50u },
   { 732515297u, 554264481u, 116u },
   { 732515297u, 554264481u, 116u },
   { 732515297u, 554264481u, 116u },
   { 2681625681u, 3960593856u, 739u },
   { 2546105329u, 1679561907u, 1875u },
   { 533319201u, 2248012685u, 5393u },
   { 533319201u, 2248012685u, 5393u },
   { 996692113u, 2949507147u, 16011u },
   { 996692113u, 2949507147u, 16011u },
   { 2616068761u, 328479117u, 198156u },
   { 1411295841u, 761797252u, 229581u }
};
#else
mp_limb_t flint_fmpz_pseudosquares[][2] =
{
   { 17, 0 },
   { 73, 0 },
   { 241, 0 },
   { 1009, 0 },
   { 2641, 0 },
   { 8089, 0 },
   { 18001, 0 },
   { 53881, 0 },
   { 87481, 0 },
   { 117049, 0 },
   { 515761, 0 },
   { 1083289, 0 },
   { 3206641, 0 },
   { 3818929, 0 },
   { 9257329, 0 },
   { 22000801, 0 },
   { 48473881, 0 },
   { 48473881, 0 },
   { 175244281, 0 },
   { 427733329, 0 },
   { 427733329, 0 },
   { 898716289u, 0 },
   { 2805544681u, 0 },
   { 2805544681u, 0 },
   { 2805544681u, 0 },
   { 10310263441u, 0 },
   { 23616331489u, 0 },
   { 85157610409u, 0 },
   { 85157610409u, 0 },
   { 196265095009u, 0 },
   { 196265095009u, 0 },
   { 2871842842801u, 0 },
   { 2871842842801u, 0 },
   { 2871842842801u, 0 },
   { 26250887023729u, 0 },
   { 26250887023729u, 0 },
   { 112434732901969u, 0 },
   { 112434732901969u, 0 },
   { 112434732901969u, 0 },
   { 178936222537081u, 0 },
   { 178936222537081u, 0 },
   { 696161110209049u, 0 },
   { 696161110209049u, 0 },
   { 2854909648103881u, 0 },
   { 6450045516630769u, 0 },
   { 6450045516630769u, 0 },
   { 11641399247947921u, 0 },
   { 11641399247947921u, 0 },
   { 190621428905186449u, 0 },
   { 196640248121928601u, 0 },
   { 712624335095093521u, 0 },
   { 1773855791877850321u, 0 },
   { 2327687064124474441u, 0 },
   { 6384991873059836689u, 0 },
   { 8019204661305419761u, 0 },
   { 10198100582046287689u, 0 },
   { 10198100582046287689u, 0 },
   { 10198100582046287689u, 0 },
   { 14508056099771532121u, 3u },
   { 6022180988239908185u, 11u },
   { 17043829275960758433u, 28u },
   { 14326875581237116289u, 50u },
   { 14326875581237116289u, 50u },
   { 2380547819961928673u, 116u },
   { 2380547819961928673u, 116u },
   { 2380547819961928673u, 116u },
   { 17010621086940159057u, 739u },
   { 7213663464718498801u, 1875u },
   { 9655140963601468961u, 5393u },
   { 9655140963601468961u, 5393u },
   { 12668036736679956625u, 16011u },
   { 12668036736679956625u, 16011u },
   { 1410807067550026393u, 198156u },
   { 3271894284933966433u, 229581u }
};
#endif

#define FLINT_NUM_FMPZ_PSEUDOSQUARES 74

void fmpz_set_pseudosquare(fmpz_t f, unsigned int i)
{
#ifndef FLINT64
   if (i < 25)
      fmpz_set_ui(f, flint_fmpz_pseudosquares[i][0]);
   else if (i < 58)
   {
      fmpz_set_ui(f, flint_fmpz_pseudosquares[i][1]);
      fmpz_mul_2exp(f, f, 32);
      fmpz_add_ui(f, f, flint_fmpz_pseudosquares[i][0]);
   } else if (i < FLINT_NUM_FMPZ_PSEUDOSQUARES)
   {
      fmpz_set_ui(f, flint_fmpz_pseudosquares[i][2]);
      fmpz_mul_2exp(f, f, 32);
      fmpz_add_ui(f, f, flint_fmpz_pseudosquares[i][1]);
      fmpz_mul_2exp(f, f, 32);
      fmpz_add_ui(f, f, flint_fmpz_pseudosquares[i][1]);
   }
#else
   if (i < 58)
      fmpz_set_ui(f, flint_fmpz_pseudosquares[i][0]);
   else if (i < FLINT_NUM_FMPZ_PSEUDOSQUARES)
   {
      fmpz_set_ui(f, flint_fmpz_pseudosquares[i][1]);
      fmpz_mul_2exp(f, f, 64);
      fmpz_add_ui(f, f, flint_fmpz_pseudosquares[i][0]);
   }
#endif
   else
   {
      flint_printf("Exception (fmpz_set_pseudosquare). Index too large.\n");
      flint_abort();
   }
}

int fmpz_is_prime_pseudosquare(const fmpz_t n)
{
    unsigned int i, j, m1;
    mp_limb_t p, B, mod8;
    fmpz_t NB, f, exp, mod, nm1;
    int ret;
    const mp_limb_t * primes;

    ret = -1; /* silence compiler warning (not set when aborting) */

    if (fmpz_sgn(n) <= 0)
       return 0;

    if (fmpz_size(n) == 1)
       return n_is_prime_pseudosquare(fmpz_get_ui(n));

    primes = n_primes_arr_readonly(FLINT_PSEUDOSQUARES_CUTOFF + 1);

    for (i = 0; i < FLINT_PSEUDOSQUARES_CUTOFF; i++)
    {
        p = primes[i];
        if (fmpz_fdiv_ui(n, p) == 0)
           return 0;
    }

    fmpz_init(NB);
    fmpz_init(f);
    fmpz_init(exp);
    fmpz_init(mod);
    fmpz_init(nm1);

    B  = primes[FLINT_PSEUDOSQUARES_CUTOFF];
    fmpz_sub_ui(nm1, n, 1);
    fmpz_fdiv_q_ui(NB, nm1, B);
    fmpz_add_ui(NB, NB, 1);

    m1 = 0;

    for (i = 0; i < FLINT_NUM_FMPZ_PSEUDOSQUARES; i++)
    {
       fmpz_set_pseudosquare(f, i);
       if (fmpz_cmp(f, NB) > 0)
          break;
    }

    if (i == FLINT_NUM_FMPZ_PSEUDOSQUARES)
    {
       ret = -1;
       goto cleanup;
    }

    fmpz_fdiv_q_2exp(exp, nm1, 1);

    for (j = 0; j <= i; j++)
    {
        fmpz_set_ui(mod, primes[j]);
        fmpz_powm(mod, mod, exp, n);
        if (!fmpz_is_one(mod) && fmpz_cmp(mod, nm1) != 0)
        {
           ret = 0;
           goto cleanup;
        }
        if (fmpz_cmp(mod, nm1) == 0)
           m1 = 1;
    }

    mod8 = fmpz_fdiv_ui(n, 8);

    if ((mod8 == 3) || (mod8 == 7))
    {
       ret = 1;
       goto cleanup;
    }

    if (mod8 == 5)
    {
        fmpz_set_ui(mod, 2);
        fmpz_powm(mod, mod, exp, n);
        if (fmpz_cmp(mod, nm1) == 0)
        {
           ret = 1;
           goto cleanup;
        }
        flint_printf("Whoah, ");
        fmpz_print(n);
        flint_printf("is a probable prime, but not prime, please report!!\n");
        flint_abort();
    }
    else
    {
        if (m1)
        {
            ret = 1;
            goto cleanup;
        }

        for (j = i + 1; j < FLINT_NUM_FMPZ_PSEUDOSQUARES + 1; j++)
        {
            fmpz_set_ui(mod, primes[j]);
            fmpz_powm(mod, mod, exp, n);
            if (fmpz_cmp(mod, nm1) == 0)
            {
               ret = 1;
               goto cleanup;
            }
            if (!fmpz_is_one(mod))
            {
                flint_printf("Whoah, ");
                fmpz_print(n);
                flint_printf("is a probable prime, but not prime, please report!!\n");
                flint_abort();
            }
        }
        flint_printf("Whoah, ");
        fmpz_print(n);
        flint_printf("is a probable prime, but not prime, please report!!\n");
        flint_abort();
    }

cleanup:

    fmpz_clear(NB);
    fmpz_clear(f);
    fmpz_clear(exp);
    fmpz_clear(mod);
    fmpz_clear(nm1);

    return ret;
}
