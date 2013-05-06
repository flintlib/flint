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
#include <stdio.h>
#define ulong unsigned long
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

#define n_pp1_set(x1, y1, x2, y2) \
   do {                           \
      x1 = x2;                    \
      y1 = y2;                    \
   } while (0)

#define n_pp1_set_ui(x, norm, c) \
   do {                          \
      x = (c << norm);           \
   } while (0)

void n_pp1_print(mp_limb_t x, mp_limb_t y, ulong norm)
{
   if (norm)
   {
      x >>= norm;
      y >>= norm;
   } 

   printf("[%lu, %lu]", x, y);
}

#define n_pp1_2k(x, y, n, ninv, x0, norm)       \
   do {                                         \
      const mp_limb_t two = (2UL << norm);      \
      y = n_mulmod_preinv(y, x, n, ninv, norm); \
      y = n_submod(y, x0, n);                   \
      x = n_mulmod_preinv(x, x, n, ninv, norm); \
      x = n_submod(x, two, n);                  \
   } while (0)

#define n_pp1_2kp1(x, y, n, ninv, x0, norm)     \
   do {                                         \
      const mp_limb_t two = (2UL << norm);      \
      x = n_mulmod_preinv(x, y, n, ninv, norm); \
      x = n_submod(x, x0, n);                   \
      y = n_mulmod_preinv(y, y, n, ninv, norm); \
      y = n_submod(y, two, n);                  \
   } while (0)

void n_pp1_pow_ui(mp_limb_t * x, mp_limb_t * y, ulong exp, 
                    mp_limb_t n, mp_limb_t ninv, ulong norm)
{
   const mp_limb_t x0 = *x;
   const mp_limb_t two = (2UL << norm);
   ulong bit = ((1UL << FLINT_BIT_COUNT(exp)) >> 2);

   (*y) = n_mulmod_preinv(*x, *x, n, ninv, norm);
   (*y) = n_submod(*y, two, n);
   
   while (bit)
   {
      if (exp & bit)
         n_pp1_2kp1(*x, *y, n, ninv, x0, norm);
      else
         n_pp1_2k(*x, *y, n, ninv, x0, norm);

      bit >>= 1;
   }
}

mp_limb_t n_pp1_factor(mp_limb_t n, mp_limb_t x, ulong norm)
{
   if (norm)
   {
      n >>= norm;
      x >>= norm;
   }
   
   x = n_submod(x, 2, n);
   if (x == 0)
      return 0;

   return n_gcd(n, x);
}

mp_limb_t n_pp1_find_power(mp_limb_t * x, mp_limb_t * y, 
                  ulong p, mp_limb_t n, mp_limb_t ninv, ulong norm)
{
   mp_limb_t factor;
   
   do
   {
      n_pp1_pow_ui(x, y, p, n, ninv, norm);
      factor = n_pp1_factor(n, *x, norm);
   } while (factor == 1);

   return factor;
}

mp_limb_t n_factor_pp1(mp_limb_t n, ulong B1, ulong c)
{
   long i, j;
   mp_limb_t factor = 0;
   mp_limb_t x, y, oldx, oldy, ninv;
   ulong pr, oldpr, sqrt, bits0, norm;
   n_primes_t iter;

   if ((n % 2) == 0)
      return 2;

   n_primes_init(iter);

   sqrt = n_sqrt(B1);
   bits0 = FLINT_BIT_COUNT(B1);

   count_leading_zeros(norm, n);
   n <<= norm;

   ninv = n_preinvert_limb(n);
   
   n_pp1_set_ui(x, norm, c);
   
   /* mul by various prime powers */   
   pr = 0;
   oldpr = 0;

   for (i = 0; pr < B1; )
   {
      j = i + 1024;
      oldpr = pr;
      n_pp1_set(oldx, oldy, x, y);
         
      for ( ; i < j; i++)
      {
         pr = n_primes_next(iter);
         if (pr < sqrt)
         {
            ulong bits = FLINT_BIT_COUNT(pr);
            ulong exp = bits0 / bits;
            n_pp1_pow_ui(&x, &y, n_pow(pr, exp), n, ninv, norm);
         } else
            n_pp1_pow_ui(&x, &y, pr, n, ninv, norm);
      }
      
      factor = n_pp1_factor(n, x, norm);
      if (factor == 0)
         break;
      if (factor != 1)
         goto cleanup;
   }

   if (pr < B1) /* factor = 0 */
   {
      n_primes_jump_after(iter, oldpr);
      n_pp1_set(x, y, oldx, oldy);  

      do
      {
         pr = n_primes_next(iter);
         n_pp1_set(oldx, oldy, x, y);
         if (pr < sqrt)
         {
            ulong bits = FLINT_BIT_COUNT(pr);
            ulong exp = bits0 / bits;
            n_pp1_pow_ui(&x, &y, n_pow(pr, exp), n, ninv, norm);
         } else
            n_pp1_pow_ui(&x, &y, pr, n, ninv, norm);

         factor = n_pp1_factor(n, x, norm);
         if (factor == 0)
            break;
         if (factor != 1)
            goto cleanup;
      } while (1);
   } else
      goto cleanup;

   /* factor still 0 */
   factor = n_pp1_find_power(&oldx, &oldy, pr, n, ninv, norm);

cleanup:

   n_primes_clear(iter);

   return factor;
}
