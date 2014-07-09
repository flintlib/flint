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
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

int fmpz_is_prime_lenstra3(fmpz_t F, fmpz * r, const fmpz_t n, 
                                                     mp_ptr pk1, slong num_pk1)
{
   slong i;
   flint_rand_t state;
   fmpz_mod_poly_t f, g, h, h1, h2, h3, finv;
   const slong k = 3; /* n^3 - 1 */
   fmpz_mod_poly_frobenius_powers_t pow;
   fmpz_t exp, exp2, Ft;
   int res = 1;

   /* multiply primes */
   fmpz_init(Ft);
   fmpz_set_ui(F, 1);
   fmpz_set_ui(Ft, 1);

   for (i = 0; i < num_pk1; i++)
      fmpz_mul_ui(Ft, Ft, pk1[i]);
   
   /* find monic irreducible poly f of degree k in Z/nZ[x] */
   flint_randinit(state);

   fmpz_mod_poly_init2(f, n, k + 1);
   fmpz_mod_poly_set_coeff_ui(f, k, 1);

   do {
      for (i = 0; i < k; i++)
         fmpz_randm(f->coeffs + i, state, n);
   } while (!fmpz_mod_poly_is_irreducible(f));

   /* compute random poly g of length k */
   fmpz_mod_poly_init2(g, n, k);

   for (i = 0; i < k; i++)
      fmpz_randm(g->coeffs + i, state, n);
   
   _fmpz_mod_poly_set_length(g, k);

   /* precompute frobenius powers for use up to x^(n^k) */
   fmpz_mod_poly_init(finv, n);

   fmpz_mod_poly_reverse(finv, f, f->length);
   fmpz_mod_poly_inv_series_newton(finv, finv, f->length);

   fmpz_mod_poly_frobenius_powers_precomp(pow, f, finv, k);

   /* test g^(n^k - 1) = 1 mod f */
   fmpz_mod_poly_init(h, n);

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

      We compute g^(n - 1) mod f and g^((n^2 + n + 1)/ q) and compose
   */
   fmpz_mod_poly_init(h1, n);
   
   /* compute x^(n - 1) */
   fmpz_mod_poly_frobenius_power(h, pow, f, 1);
   fmpz_mod_poly_set_coeff_ui(h1, 1, 1);
   fmpz_mod_poly_set_coeff_ui(h1, 0, 0);
   _fmpz_mod_poly_set_length(h1, 2);
   fmpz_mod_poly_invmod(h1, h1, f);
   fmpz_mod_poly_mulmod_preinv(h, h, h1, f, finv);
   
   /* compute x^((n^2 + n + 1) / q) for each prime in pk1 */
   fmpz_init(exp);
   fmpz_init(exp2);
   
   fmpz_mod_poly_init(h2, n);
   fmpz_mod_poly_init(h3, n);

   /* compute x^((n^2 + n + 1)/Ft) */
   fmpz_mul(exp, n, n);
   fmpz_add(exp, exp, n);
   fmpz_add_ui(exp, exp, 1);
   fmpz_tdiv_q(exp, exp, Ft);
   fmpz_mod_poly_powmod_x_fmpz_preinv(h3, exp2, f, finv);
       
   /* raise to powers F/q and compose with g to get g^((n^2 + n + 1)/q) */
   for (i = 0; i < num_pk1; i++)
   {
      fmpz_divexact_ui(exp2, Ft, pk1[i]);

      fmpz_mod_poly_powmod_x_fmpz_preinv(h2, exp2, f, finv);
      fmpz_mod_poly_compose_mod(h2, h2, h3, f);
      fmpz_mod_poly_compose_mod(h1, h, h2, f);
      fmpz_mod_poly_compose_mod(h1, h1, g, f);

      if (!fmpz_mod_poly_is_one(h1))
         fmpz_mul_ui(F, F, pk1[i]);
   }

   fmpz_mod_poly_clear(h1);
   fmpz_mod_poly_clear(h2);
   fmpz_mod_poly_clear(h3);

   fmpz_clear(exp);
   fmpz_clear(exp2);

   /* set residues to n^j mod F for j = 0, k - 1 */
   fmpz_set_ui(r + 0, 1);
   fmpz_mod(r + 1, n, F);
   fmpz_mul(r + 2, r + 1, r + 1);
   fmpz_mod(r + 2, r + 2, F);

cleanup:

   fmpz_mod_poly_frobenius_powers_clear(pow);

   fmpz_mod_poly_clear(finv);
   fmpz_mod_poly_clear(f);
   fmpz_mod_poly_clear(g);
   fmpz_mod_poly_clear(h);

   fmpz_clear(Ft);

   flint_randclear(state);

   return res;
}

int fmpz_is_prime_lenstra4(fmpz_t F, fmpz * r, const fmpz_t n, 
                                                     mp_ptr pk1, slong num_pk1)
{
   slong i;
   flint_rand_t state;
   fmpz_mod_poly_t f, g, h, h1, h2, h3, finv;
   const slong k = 4; /* n^4 - 1 */
   fmpz_mod_poly_frobenius_powers_t pow;
   fmpz_t exp, exp2, Ft;
   int res = 1;

   /* multiply primes */
   fmpz_init(Ft);
   fmpz_set_ui(F, 1);
   fmpz_set_ui(Ft, 1);

   for (i = 0; i < num_pk1; i++)
      fmpz_mul_ui(Ft, Ft, pk1[i]);
   
   /* find monic irreducible poly f of degree k in Z/nZ[x] */
   flint_randinit(state);

   fmpz_mod_poly_init2(f, n, k + 1);
   fmpz_mod_poly_set_coeff_ui(f, k, 1);

   do {
      for (i = 0; i < k; i++)
         fmpz_randm(f->coeffs + i, state, n);
   } while (!fmpz_mod_poly_is_irreducible(f));

   /* compute random poly g of length k */
   fmpz_mod_poly_init2(g, n, k);

   for (i = 0; i < k; i++)
      fmpz_randm(g->coeffs + i, state, n);
   
   _fmpz_mod_poly_set_length(g, k);

   /* precompute frobenius powers for use up to x^(n^k) */
   fmpz_mod_poly_init(finv, n);

   fmpz_mod_poly_reverse(finv, f, f->length);
   fmpz_mod_poly_inv_series_newton(finv, finv, f->length);

   fmpz_mod_poly_frobenius_powers_precomp(pow, f, finv, k);

   /* test g^(n^k - 1) = 1 mod f */
   fmpz_mod_poly_init(h, n);

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

      We compute g^(n^2 - 1) mod f and g^((n^2 + 1)/ q) and compose
   */
   fmpz_mod_poly_init(h1, n);
   
   /* compute x^(n^2 - 1) */
   fmpz_mod_poly_frobenius_power(h, pow, f, 2);
   fmpz_mod_poly_set_coeff_ui(h1, 1, 1);
   fmpz_mod_poly_set_coeff_ui(h1, 0, 0);
   _fmpz_mod_poly_set_length(h1, 2);
   fmpz_mod_poly_invmod(h1, h1, f);
   fmpz_mod_poly_mulmod_preinv(h, h, h1, f, finv);
   
   /* compute x^((n^2 + 1) / q) for each prime in pk1 */
   fmpz_init(exp);
   fmpz_init(exp2);
   
   fmpz_mod_poly_init(h2, n);
   fmpz_mod_poly_init(h3, n);

   /* compute x^((n^2 + 1)/Ft) */
   fmpz_mul(exp, n, n);
   fmpz_add_ui(exp, exp, 1);
   fmpz_tdiv_q(exp, exp, Ft);
   fmpz_mod_poly_powmod_x_fmpz_preinv(h3, exp2, f, finv);
       
   /* raise to powers F/q and compose with g to get g^((n^2 + 1)/q) */
   for (i = 0; i < num_pk1; i++)
   {
      fmpz_divexact_ui(exp2, Ft, pk1[i]);

      fmpz_mod_poly_powmod_x_fmpz_preinv(h2, exp2, f, finv);
      fmpz_mod_poly_compose_mod(h2, h2, h3, f);
      fmpz_mod_poly_compose_mod(h1, h, h2, f);
      fmpz_mod_poly_compose_mod(h1, h1, g, f);

      if (!fmpz_mod_poly_is_one(h1))
         fmpz_mul_ui(F, F, pk1[i]);
   }

   fmpz_mod_poly_clear(h1);
   fmpz_mod_poly_clear(h2);
   fmpz_mod_poly_clear(h3);

   fmpz_clear(exp);
   fmpz_clear(exp2);

   /* set residues to n^j mod F for j = 0, k - 1 */
   fmpz_set_ui(r + 0, 1);
   fmpz_mod(r + 1, n, F);
   fmpz_mul(r + 2, r + 1, r + 1);
   fmpz_mod(r + 2, r + 2, F);
   fmpz_mul(r + 3, r + 2, r + 2);
   fmpz_mod(r + 3, r + 3, F);

cleanup:

   fmpz_mod_poly_frobenius_powers_clear(pow);

   fmpz_mod_poly_clear(finv);
   fmpz_mod_poly_clear(f);
   fmpz_mod_poly_clear(g);
   fmpz_mod_poly_clear(h);

   fmpz_clear(Ft);

   flint_randclear(state);

   return res;
}

int fmpz_is_prime_lenstra6(fmpz_t F, fmpz * r, const fmpz_t n, 
                                                     mp_ptr pk1, slong num_pk1)
{
   slong i;
   flint_rand_t state;
   fmpz_mod_poly_t f, g, h, h1, h2, h3, finv;
   const slong k = 6; /* n^6 - 1 */
   fmpz_mod_poly_frobenius_powers_t pow;
   fmpz_t exp, exp2, Ft;
   int res = 1;

   /* multiply primes */
   fmpz_init(Ft);
   fmpz_set_ui(F, 1);
   fmpz_set_ui(Ft, 1);

   for (i = 0; i < num_pk1; i++)
      fmpz_mul_ui(Ft, Ft, pk1[i]);
   
   /* find monic irreducible poly f of degree k in Z/nZ[x] */
   flint_randinit(state);

   fmpz_mod_poly_init2(f, n, k + 1);
   fmpz_mod_poly_set_coeff_ui(f, k, 1);

   do {
      for (i = 0; i < k; i++)
         fmpz_randm(f->coeffs + i, state, n);
   } while (!fmpz_mod_poly_is_irreducible(f));

   /* compute random poly g of length k */
   fmpz_mod_poly_init2(g, n, k);

   for (i = 0; i < k; i++)
      fmpz_randm(g->coeffs + i, state, n);
   
   _fmpz_mod_poly_set_length(g, k);

   /* precompute frobenius powers for use up to x^(n^k) */
   fmpz_mod_poly_init(finv, n);

   fmpz_mod_poly_reverse(finv, f, f->length);
   fmpz_mod_poly_inv_series_newton(finv, finv, f->length);

   fmpz_mod_poly_frobenius_powers_precomp(pow, f, finv, k);

   /* test g^(n^k - 1) = 1 mod f */
   fmpz_mod_poly_init(h, n);

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

      We compute g^(n - 1) mod f and g^((n^2 - n + 1)/ q) and compose
   */
   fmpz_mod_poly_init(h1, n);
   fmpz_mod_poly_init(h2, n);
   
   /* compute x^(n^3 - 1) */
   fmpz_mod_poly_frobenius_power(h, pow, f, 3);
   fmpz_mod_poly_set_coeff_ui(h1, 1, 1);
   fmpz_mod_poly_set_coeff_ui(h1, 0, 0);
   _fmpz_mod_poly_set_length(h1, 2);
   fmpz_mod_poly_invmod(h1, h1, f);
   fmpz_mod_poly_mulmod_preinv(h, h, h1, f, finv);
   
   /* compute x^(n + 1) */
   fmpz_mod_poly_frobenius_power(h2, pow, f, 1);
   fmpz_mod_poly_set_coeff_ui(h1, 1, 1);
   fmpz_mod_poly_set_coeff_ui(h1, 0, 0);
   _fmpz_mod_poly_set_length(h1, 2);
   fmpz_mod_poly_mulmod_preinv(h2, h2, h1, f, finv);
   
   /* compute x^((n^3 - 1)(n + 1)) */
   fmpz_mod_poly_compose_mod(h, h, h2, f);

   /* compute x^((n^2 - n + 1) / q) for each prime in pk1 */
   fmpz_init(exp);
   fmpz_init(exp2);
   
   fmpz_mod_poly_init(h3, n);

   /* compute x^((n^2 - n + 1)/Ft) */
   fmpz_mul(exp, n, n);
   fmpz_sub(exp, exp, n);
   fmpz_add_ui(exp, exp, 1);
   fmpz_tdiv_q(exp, exp, Ft);
   fmpz_mod_poly_powmod_x_fmpz_preinv(h3, exp2, f, finv);
       
   /* raise to powers F/q and compose with g to get g^((n^2 - n + 1)/q) */
   for (i = 0; i < num_pk1; i++)
   {
      fmpz_divexact_ui(exp2, Ft, pk1[i]);

      fmpz_mod_poly_powmod_x_fmpz_preinv(h2, exp2, f, finv);
      fmpz_mod_poly_compose_mod(h2, h2, h3, f);
      fmpz_mod_poly_compose_mod(h1, h, h2, f);
      fmpz_mod_poly_compose_mod(h1, h1, g, f);

      if (!fmpz_mod_poly_is_one(h1))
         fmpz_mul_ui(F, F, pk1[i]);
   }

   fmpz_mod_poly_clear(h1);
   fmpz_mod_poly_clear(h2);
   fmpz_mod_poly_clear(h3);

   fmpz_clear(exp);
   fmpz_clear(exp2);

   /* set residues to n^j mod F for j = 0, k - 1 */
   fmpz_set_ui(r + 0, 1);
   fmpz_mod(r + 1, n, F);
   for (i = 2; i < 6; i++)
   {
      fmpz_mul(r + i, r + i - 1, r + i - 1);
      fmpz_mod(r + i, r + i, F);
   }

cleanup:

   fmpz_mod_poly_frobenius_powers_clear(pow);

   fmpz_mod_poly_clear(finv);
   fmpz_mod_poly_clear(f);
   fmpz_mod_poly_clear(g);
   fmpz_mod_poly_clear(h);

   fmpz_clear(Ft);

   flint_randclear(state);

   return res;
}