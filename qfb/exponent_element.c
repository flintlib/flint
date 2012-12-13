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

    Copyright (C) 2012 William Hart

******************************************************************************/

#undef ulong /* prevent clash with stdlib */
#include <stdlib.h>
#define ulong unsigned long
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qfb.h"

/*
   find which power of the base is the exponent of f
*/
ulong find_power(qfb_t f, fmpz_t n, ulong base)
{
   ulong s = 1;
      
   do
   {
      qfb_pow_ui(f, f, n, base);
      s *= base;
   } while (!qfb_is_principal_form(f, n));

   return s;
}

int qfb_exponent_element(fmpz_t exponent, qfb_t f, fmpz_t n, long iters)
{
   long i, j;
   qfb_t pow, oldpow;
   ulong pr, oldpr;
   
   if (qfb_is_principal_form(f, n))
   {
      fmpz_set_ui(exponent, 1);
      return 1;
   }

   qfb_init(pow);
   qfb_init(oldpow);
   
   qfb_set(pow, f);
   
   /* raise to various powers of small primes */
   qfb_set(oldpow, pow);
   qfb_pow_ui(pow, oldpow, n, 4096); /* 2^12 */
   if (qfb_is_principal_form(pow, n))
   {
      ulong s = find_power(oldpow, n, 2);
      fmpz_set_ui(exponent, s);

      return 1;
   }
   
   qfb_set(oldpow, pow);
   qfb_pow_ui(pow, oldpow, n, 59049); /* 3^10 */
   if (qfb_is_principal_form(pow, n))
   {
      ulong s = find_power(oldpow, n, 3);
      qfb_pow_ui(pow, f, n, s);
      qfb_exponent_element(exponent, pow, n, iters);
      fmpz_mul_ui(exponent, exponent, s);
      
      return 1;
   }
   
   qfb_set(oldpow, pow);
   qfb_pow_ui(pow, oldpow, n, 390625); /* 5^8 */
   if (qfb_is_principal_form(pow, n))
   {
      ulong s = find_power(oldpow, n, 5);
      qfb_pow_ui(pow, f, n, s);
      qfb_exponent_element(exponent, pow, n, iters);
      fmpz_mul_ui(exponent, exponent, s);
      
      return 1;
   }

   qfb_set(oldpow, pow);
   qfb_pow_ui(pow, oldpow, n, 117649); /* 7^6 */
   if (qfb_is_principal_form(pow, n))
   {
      ulong s = find_power(oldpow, n, 7);
      qfb_pow_ui(pow, f, n, s);
      qfb_exponent_element(exponent, pow, n, iters);
      fmpz_mul_ui(exponent, exponent, s);
      
      return 1;
   }

   qfb_set(oldpow, pow);
   qfb_pow_ui(pow, oldpow, n, 14641); /* 11^4 */
   if (qfb_is_principal_form(pow, n))
   {
      ulong s = find_power(oldpow, n, 11);
      qfb_pow_ui(pow, f, n, s);
      qfb_exponent_element(exponent, pow, n, iters);
      fmpz_mul_ui(exponent, exponent, s);
      
      return 1;
   }

   qfb_set(oldpow, pow);
   qfb_pow_ui(pow, oldpow, n, 169); /* 13^2 */
   if (qfb_is_principal_form(pow, n))
   {
      ulong s = find_power(oldpow, n, 13);
      qfb_pow_ui(pow, f, n, s);
      qfb_exponent_element(exponent, pow, n, iters);
      fmpz_mul_ui(exponent, exponent, s);
      
      return 1;
   }

   /* raise to prime powers (squared) until the identity is found */
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

      /* identity is found, compute exponent recursively */
      if (qfb_is_principal_form(pow, n))
      {
         qfb_set(pow, oldpow);
         pr = oldpr;
         while (1)
         {
            pr = n_nextprime(pr, 0);
            
            qfb_pow_ui(pow, pow, n, pr);
            if (qfb_is_principal_form(pow, n))
            {
               qfb_pow_ui(pow, f, n, pr);
               qfb_exponent_element(exponent, pow, n, iters);
               fmpz_mul_ui(exponent, exponent, pr);
               return 1;
            }

            qfb_pow_ui(pow, pow, n, pr);
            if (qfb_is_principal_form(pow, n))
            {
               qfb_pow_ui(pow, f, n, pr*pr);
               qfb_exponent_element(exponent, pow, n, iters);
               fmpz_mul_ui(exponent, exponent, pr*pr);
               return 1;
            }
         }
      }
   }

   qfb_clear(pow);
   qfb_clear(oldpow);

   return 0;
}
