/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

/* bits of n, B1, count */
slong n_factor_pp1_table[][2] = {
    {2784, 5}, {1208, 2}, {2924, 3},
    { 286, 5}, {  58, 5}, {  61, 4}, { 815, 2},
    { 944, 2}, {  61, 3}, {   0, 0}, {   0, 0},
    {   0, 0}, {   0, 0}, {   0, 0}, {   0, 0},
    {   0, 0}, {   0, 0}, {   0, 0}, {   0, 0},
    { 606, 1}, {2403, 1}, {2524, 1}, {2924, 1},
    {3735, 2}, { 669, 2}, {6092, 3}, {2179, 3},
    {3922, 3}, {6717, 4}, {4119, 4}, {2288, 4},
    {9004, 3}, {9004, 3}, {9004, 3}};

#define n_pp1_set(x1, y1, x2, y2) \
   do {                           \
      x1 = x2;                    \
      y1 = y2;                    \
   } while (0)

#define n_pp1_set_ui(x, norm, c) \
   do {                          \
      x = (c << norm);           \
   } while (0)

#if 0
/* For debugging */
void n_pp1_print(ulong x, ulong y, ulong norm)
{
   if (norm)
   {
      x >>= norm;
      y >>= norm;
   }

   flint_printf("[%wu, %wu]", x, y);
}
#endif

#define n_pp1_2k(x, y, n, ninv, x0, norm)       \
   do {                                         \
      const ulong xxtwo = (UWORD(2) << norm);   \
      y = n_mulmod_preinv(y, x, n, ninv, norm); \
      y = n_submod(y, x0, n);                   \
      x = n_mulmod_preinv(x, x, n, ninv, norm); \
      x = n_submod(x, xxtwo, n);                \
   } while (0)

#define n_pp1_2kp1(x, y, n, ninv, x0, norm)     \
   do {                                         \
      const ulong xxtwo = (UWORD(2) << norm);   \
      x = n_mulmod_preinv(x, y, n, ninv, norm); \
      x = n_submod(x, x0, n);                   \
      y = n_mulmod_preinv(y, y, n, ninv, norm); \
      y = n_submod(y, xxtwo, n);                \
   } while (0)

static void n_pp1_pow_ui(ulong * x, ulong * y, ulong exp,
                    ulong n, ulong ninv, ulong norm)
{
   const ulong x0 = *x;
   const ulong two = (UWORD(2) << norm);
   ulong bit = ((UWORD(1) << FLINT_BIT_COUNT(exp)) >> 2);

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

static ulong n_pp1_factor(ulong n, ulong x, ulong norm)
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

static ulong n_pp1_find_power(ulong * x, ulong * y,
                  ulong p, ulong n, ulong ninv, ulong norm)
{
   ulong factor;

   do
   {
      n_pp1_pow_ui(x, y, p, n, ninv, norm);
      factor = n_pp1_factor(n, *x, norm);
   } while (factor == 1);

   return factor;
}

ulong n_factor_pp1(ulong n, ulong B1, ulong c)
{
   slong i, j;
   ulong factor = 0;
   ulong x, y = 0, oldx, oldy, ninv;
   ulong pr, oldpr, sqrt, bits0, norm;
   n_primes_t iter;

   if ((n % 2) == 0)
      return 2;

   n_primes_init(iter);

   sqrt = n_sqrt(B1);
   bits0 = FLINT_BIT_COUNT(B1);

   norm = flint_clz(n);
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
   {
      factor = 0;
      goto cleanup;
   }

   /* factor still 0 */
   factor = n_pp1_find_power(&oldx, &oldy, pr, n, ninv, norm);

cleanup:

   n_primes_clear(iter);

   return factor;
}

ulong n_factor_pp1_wrapper(ulong n)
{
   slong bits = FLINT_BIT_COUNT(n);
   ulong B1;
   slong count, i;
   flint_rand_t state;

   /* silently fail if trial factoring would always succeed */
   if (bits < 31)
       return 0;

   B1 = n_factor_pp1_table[bits - 31][0];
   count = n_factor_pp1_table[bits - 31][1];

   flint_rand_init(state);

   for (i = 0; i < count; i++)
   {
       ulong factor;
       factor = n_factor_pp1(n, B1, n_randint(state, n - 3) + 3);
       if (factor != 0)
       {
           flint_rand_clear(state);
           return factor;
       }
   }

   flint_rand_clear(state);
   return 0;
}
