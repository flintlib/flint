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
#include <time.h>
#include <stdlib.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

int fmpz_is_prime_lenstra(fmpz_t F, fmpz * r, const fmpz_t n, 
                                                     mp_ptr pk1, slong num_pk1, slong k)
{
   slong i, bits = fmpz_bits(n);
   flint_rand_t state;
   fmpz_mod_poly_t f, g, h, h1, h2, h3, finv;
   fmpz_mod_poly_frobenius_powers_2exp_t pow;
   fmpz_t exp, exp2, Ft, a1, a2, b, p, R, nfac, base;
   fmpz_factor_t fac;
   int res = 1, fac_found, d;

   fmpz_set_ui(F, 1);
   
   if (num_pk1 == 0)
      goto end;
   
   fmpz_factor_init(fac);

   /* compute polynomial */
   fmpz_init(exp);
   fmpz_init(exp2);
   fmpz_init(R);
   fmpz_init(p);
   fmpz_init(nfac);

   switch (k)
   {
   case 3:
      /* n^2 + n + 1 */
      fmpz_mul(exp, n, n);
      fmpz_add(exp, exp, n);
      break;
   case 4:
      /* n^2 + 1 */
      fmpz_mul(exp, n, n);
      break;
   case 8:
      /* n^4 + 1 */
      fmpz_mul(exp, n, n);
      fmpz_mul(exp, exp, exp);
      break;
   case 5:
      /* n^4 + n^3 + n^2 + n + 1 */
      fmpz_mul(exp2, n, n);
      fmpz_mul(exp, exp2, exp2);
      fmpz_addmul(exp, exp2, n);
      fmpz_add(exp, exp, exp2);
      fmpz_add(exp, exp, n);
      break;
   case 10:
      /* n^4 - n^3 + n^2 - n + 1 */
      fmpz_mul(exp2, n, n);
      fmpz_mul(exp, exp2, exp2);
      fmpz_submul(exp, exp2, n);
      fmpz_add(exp, exp, exp2);
      fmpz_sub(exp, exp, n);
      break;
   case 6:
      /* n^2 - n + 1 */
      fmpz_mul(exp, n, n);
      fmpz_sub(exp, exp, n);
      break;
   case 12:
      /* n^4 - n^2 + 1 */
      fmpz_mul(exp2, n, n);
      fmpz_mul(exp, exp2, exp2);
      fmpz_sub(exp, exp, exp2);
      break;
   }

   fmpz_add_ui(exp, exp, 1);
   
   fmpz_set(R, exp);

   /* remove already found trial factors */
   for (i = 0; i < num_pk1; i++)
   {
      fmpz_set_ui(p, pk1[i]);
      d = fmpz_remove(R, R, p);
      _fmpz_factor_append(fac, p, d);
   }

   srand(time(NULL));

   /* multiply primes */
   fmpz_init(Ft);

   fmpz_set_ui(Ft, 1);

   for (i = 0; i < num_pk1; i++)
      fmpz_mul_ui(Ft, Ft, pk1[i]);

   /* try to factor further */
   do
   {
      if ((fac_found = fmpz_factor_pp1(p, R, bits*2 + 10, bits/10 + 10, rand()%100 + 3)
                    && fmpz_is_prime(p)))
      {
         d = fmpz_remove(R, R, p);
         _fmpz_factor_append(fac, p, d);
         fmpz_mul(Ft, Ft, p);
      }
   } while (fac_found);
   
   /* find monic irreducible poly f of degree k in Z/nZ[x] */
   flint_randinit(state);

   fmpz_mod_poly_init2(f, n, k + 1);
   fmpz_mod_poly_init(h, n);
   fmpz_mod_poly_init(h1, n);
   fmpz_mod_poly_init(h2, n);
   fmpz_mod_poly_init2(g, n, k);

   fmpz_set_ui(nfac, 1);
   fmpz_init_set_ui(base, 1);
      
   if ((k % 2) == 0)
   {
      /* 
         Find irreducible poly g of degree k/2 such that g(2)g(-2) is non-square in F_n 
         See "The explicit construction of irreducible polynomials over finite fields
         Stephen D. Cohen, Designs, Codes and Cryptography, 2, 169-174 (1992)
      */
      fmpz_mod_poly_set_coeff_ui(g, k/2, 1);

      fmpz_init(a1);
      fmpz_init(a2);
      fmpz_init(b);

      do {  
         do {
            /* prevent very long loop for composite n */
            fmpz_add_ui(base, base, 1);
            if (!fmpz_is_strong_probabprime(n, base))
            {
               fmpz_clear(base);
               fmpz_clear(a1);
               fmpz_clear(a2);
               fmpz_clear(b);

               res = 0;

               goto cleanup1;
            }

            for (i = 0; i < k/2; i++)
               fmpz_randm(g->coeffs + i, state, n);
            fmpz_set_ui(b, 2);
            fmpz_mod_poly_evaluate_fmpz(a1, g, b);
            fmpz_sub_ui(b, n, 2);
            fmpz_mod_poly_evaluate_fmpz(a2, g, b);
            fmpz_mul(a1, a1, a2);
            fmpz_mod(a1, a1, n);
         } while (fmpz_jacobi(a1, n) != -1 
              || (!fmpz_mod_poly_is_squarefree_f(nfac, g) && fmpz_is_one(nfac))
              || (!fmpz_mod_poly_is_irreducible_rabin_f(nfac, g) && fmpz_is_one(nfac)));

         if (!fmpz_is_one(nfac)) /* factor of n found */
            break;

         /* apply Q operator f(x) = x^k/2 g(x + 1/x) */

         fmpz_mod_poly_set_coeff_fmpz(f, k/2, g->coeffs + 0);

         fmpz_mod_poly_set_coeff_ui(h, 2, 1); /* x^2 + 1 */
         fmpz_mod_poly_set_coeff_ui(h, 0, 1);

         fmpz_mod_poly_set(h1, h);

         for (i = 1; i <= k/2; i++)
         {
            fmpz_mod_poly_shift_left(h2, h1, k/2 - i);
            fmpz_mod_poly_scalar_mul_fmpz(h2, h2, g->coeffs + i);
            fmpz_mod_poly_add(f, f, h2);
            if (i != k/2)
               fmpz_mod_poly_mul(h1, h1, h);
         }
      } while ((!fmpz_mod_poly_is_squarefree_f(nfac, f) && fmpz_is_one(nfac))
            || (!fmpz_mod_poly_is_irreducible_rabin_f(nfac, f) && fmpz_is_one(nfac)));

      fmpz_clear(a1);
      fmpz_clear(a2);
      fmpz_clear(b);
   } else /* k is odd */
   {
      fmpz_mod_poly_set_coeff_ui(f, k, 1);

      do {
         fmpz_add_ui(base, base, 1);
         if (!fmpz_is_strong_probabprime(n, base))
         {
            fmpz_clear(base);
            fmpz_clear(a1);
            fmpz_clear(a2);
            fmpz_clear(b);

            res = 0;

            goto cleanup1;
         }
            
         for (i = 0; i < k; i++)
            fmpz_randm(f->coeffs + i, state, n);
      } while (!fmpz_mod_poly_is_irreducible_rabin_f(nfac, f) && fmpz_is_one(nfac));
   }
      
   fmpz_clear(base);

   if (!fmpz_is_one(nfac)) /* factor of n found */
   {
      res = 0;

      goto cleanup1;
   }

   /* compute random poly g of length k */

   for (i = 0; i < k; i++)
      fmpz_randm(g->coeffs + i, state, n);
   
   _fmpz_mod_poly_set_length(g, k);
   _fmpz_mod_poly_normalise(g);

   /* precompute frobenius powers for use up to x^(n^k) */
   fmpz_mod_poly_init(finv, n);

   fmpz_mod_poly_reverse(finv, f, f->length);
   fmpz_mod_poly_inv_series_newton_f(nfac, finv, finv, f->length);

   if (!fmpz_is_one(nfac))
   {
      res = 0;

      goto cleanup2;
   }

   fmpz_mod_poly_frobenius_powers_2exp_precomp(pow, f, finv, k);

   /* test g^(n^k - 1) = 1 mod f */
   fmpz_mod_poly_frobenius_power(h, pow, f, k);
   fmpz_mod_poly_compose_mod(h, h, g, f);

   if (!fmpz_mod_poly_equal(h, g))
   {
      res = 0;

      goto cleanup;
   }

   /* 
      Determine usable primes.
      These are q such that g^((n^k - 1)/q) != 1 mod f 
   */

   switch (k)
   {
   case 3:
   case 4:
   case 5:
   case 8:
      /* compute x^(n - 1) if k is odd, else x^(n^2 - 1) or x^(n^4 - 1) */
      fmpz_mod_poly_frobenius_power(h, pow, f, 1 + ((k % 2) == 0) + 2*(k == 8));
      fmpz_mod_poly_set_coeff_ui(h1, 1, 1);
      fmpz_mod_poly_set_coeff_ui(h1, 0, 0);
      _fmpz_mod_poly_set_length(h1, 2);
      fmpz_mod_poly_invmod_f(nfac, h1, h1, f);
      if (!fmpz_is_one(nfac))
         goto cleanup;
      fmpz_mod_poly_mulmod_preinv(h, h, h1, f, finv);
      break;
   case 6:
   case 12:
   case 10:
      /* compute x^(n^(k/2) - 1) */
      fmpz_mod_poly_frobenius_power(h, pow, f, k/2);
      fmpz_mod_poly_set_coeff_ui(h1, 1, 1);
      fmpz_mod_poly_set_coeff_ui(h1, 0, 0);
      _fmpz_mod_poly_set_length(h1, 2);
      fmpz_mod_poly_invmod_f(nfac, h1, h1, f);
      if (!fmpz_is_one(nfac))
         goto cleanup;
      fmpz_mod_poly_mulmod_preinv(h, h, h1, f, finv);
   
      /* compute x^(n + 1) or x^(n^2 + 1) */
      fmpz_mod_poly_frobenius_power(h2, pow, f, 1 + ((k % 4) == 0));
      fmpz_mod_poly_set_coeff_ui(h1, 1, 1);
      fmpz_mod_poly_set_coeff_ui(h1, 0, 0);
      _fmpz_mod_poly_set_length(h1, 2);
      fmpz_mod_poly_mulmod_preinv(h2, h2, h1, f, finv);
   
      /* compose x^(n^(k/2) - 1) with x^(n + 1) or x^(n^2 + 1) */
      fmpz_mod_poly_compose_mod(h, h, h2, f);
      break;
   }

   /* compute x^(pol(n)/Ft) */
   fmpz_mod_poly_init(h3, n);

   fmpz_tdiv_q(exp, exp, Ft);
   fmpz_mod_poly_powmod_x_fmpz_preinv(h3, exp, f, finv);
       
   /* raise to powers Ft/q and compose with g to get g^(pol(n)/q) */
   for (i = 0; i < fac->num; i++)
   {
      fmpz_divexact(exp2, Ft, fac->p + i);

      fmpz_mod_poly_powmod_x_fmpz_preinv(h2, exp2, f, finv);
      fmpz_mod_poly_compose_mod(h2, h2, h3, f);
      fmpz_mod_poly_compose_mod(h1, h, h2, f);
      fmpz_mod_poly_compose_mod(h1, h1, g, f);

      if (!fmpz_mod_poly_is_one(h1))
      {
         fmpz_pow_ui(p, fac->p + i, fac->exp[i]);
         fmpz_mul(F, F, p);
      }
   }

   fmpz_mod_poly_clear(h3);

cleanup:
   fmpz_mod_poly_frobenius_powers_2exp_clear(pow);

cleanup2:

   fmpz_mod_poly_clear(finv);

cleanup1:

   fmpz_factor_clear(fac);

   fmpz_mod_poly_clear(f);
   fmpz_mod_poly_clear(g);
   fmpz_mod_poly_clear(h);
   fmpz_mod_poly_clear(h1);
   fmpz_mod_poly_clear(h2);

   fmpz_clear(p);
   fmpz_clear(Ft);
   fmpz_clear(exp);
   fmpz_clear(exp2);
   fmpz_clear(R);
   fmpz_clear(nfac);

   flint_randclear(state);

end:

   /* set residues to n^j mod F for j = 0, k - 1 */
   if (fmpz_is_one(F))
      fmpz_set_ui(r + 0, 0);
   else
      fmpz_set_ui(r + 0, 1);
   fmpz_mod(r + 1, n, F);
   for (i = 2; i < k; i++)
   {
      fmpz_mul(r + i, r + i - 1, r + i - 1);
      fmpz_mod(r + i, r + i, F);
   }

   return res;
}
