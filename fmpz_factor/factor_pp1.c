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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_factor.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

void pp1_set(mp_ptr x1, mp_ptr y1, 
              mp_srcptr x2, mp_srcptr y2, mp_size_t nn)
{
   mpn_copyi(x1, x2, nn);
   mpn_copyi(y1, y2, nn);
}

void pp1_set_ui(mp_ptr x, mp_size_t nn, ulong norm, ulong c)
{
   mpn_zero(x, nn);
   x[0] = (c << norm);
   if (nn > 1 && norm)
      x[1] = (c >> (FLINT_BITS - norm));
}

void pp1_print(mp_srcptr x, mp_srcptr y, mp_size_t nn, ulong norm)
{
   mp_ptr tx = flint_malloc(nn*sizeof(mp_limb_t));
   mp_ptr ty = flint_malloc(nn*sizeof(mp_limb_t));

   if (norm)
   {
      mpn_rshift(tx, x, nn, norm);
      mpn_rshift(ty, y, nn, norm);
   } else
   {
       mpn_copyi(tx, x, nn);
       mpn_copyi(ty, y, nn);
   }

   printf("["), gmp_printf("%Nd", tx, nn), printf(", "), gmp_printf("%Nd", ty, nn), printf("]");

   flint_free(tx);
   flint_free(ty);
}

void pp1_2k(mp_ptr x, mp_ptr y, mp_size_t nn, mp_srcptr n, 
            mp_srcptr ninv, mp_srcptr x0, ulong norm)
{
   flint_mpn_mulmod_preinvn(y, y, x, nn, n, ninv, norm);
   if (mpn_sub_n(y, y, x0, nn))
      mpn_add_n(y, y, n, nn);

   flint_mpn_mulmod_preinvn(x, x, x, nn, n, ninv, norm);
   if (mpn_sub_1(x, x, nn, 2UL << norm))
      mpn_add_n(x, x, n, nn);
}

void pp1_2kp1(mp_ptr x, mp_ptr y, mp_size_t nn, mp_srcptr n, 
              mp_srcptr ninv, mp_srcptr x0, ulong norm)
{
   flint_mpn_mulmod_preinvn(x, x, y, nn, n, ninv, norm);
   if (mpn_sub_n(x, x, x0, nn))
      mpn_add_n(x, x, n, nn);

   flint_mpn_mulmod_preinvn(y, y, y, nn, n, ninv, norm);
   if (mpn_sub_1(y, y, nn, 2UL << norm))
      mpn_add_n(y, y, n, nn);
}

void pp1_pow_ui(mp_ptr x, mp_ptr y, mp_size_t nn, 
                ulong exp, mp_srcptr n, mp_srcptr ninv, ulong norm)
{
   mp_limb_t t[30];
   mp_ptr x0 = t;
   ulong bit = ((1UL << FLINT_BIT_COUNT(exp)) >> 2);

   if (nn > 30)
      x0 = flint_malloc(nn*sizeof(mp_limb_t));
   mpn_copyi(x0, x, nn);

   flint_mpn_mulmod_preinvn(y, x, x, nn, n, ninv, norm);
   if (mpn_sub_1(y, y, nn, 2UL << norm))
      mpn_add_n(y, y, n, nn);

   while (bit)
   {
      if (exp & bit)
         pp1_2kp1(x, y, nn, n, ninv, x0, norm);
      else
         pp1_2k(x, y, nn, n, ninv, x0, norm);

      bit >>= 1;
   }

   if (nn > 30)
      flint_free(x0);
}

mp_size_t pp1_factor(mp_ptr factor, mp_srcptr n,
                     mp_srcptr x, mp_size_t nn, ulong norm)
{
   mp_size_t ret = 0, xn = nn;
   
   mp_ptr n2 = flint_malloc(nn*sizeof(mp_limb_t));
   mp_ptr x2 = flint_malloc(nn*sizeof(mp_limb_t));

   if (norm)
      mpn_rshift(n2, n, nn, norm);
   else
      mpn_copyi(n2, n, nn);

   if (norm)
      mpn_rshift(x2, x, nn, norm);
   else
      mpn_copyi(x2, x, nn);
   
   if (mpn_sub_1(x2, x2, nn, 2))
      mpn_add_n(x2, x2, n2, nn);

   MPN_NORM(x2, xn);

   if (xn == 0)
      goto cleanup;
   
   ret = flint_mpn_gcd_full(factor, n2, nn, x2, xn);

cleanup:

   flint_free(n2);
   flint_free(x2);

   return ret;
}

mp_size_t pp1_find_power(mp_ptr factor, mp_ptr x, mp_ptr y, mp_size_t nn, 
                          ulong p, mp_srcptr n, mp_srcptr ninv, ulong norm)
{
   mp_size_t ret;
   
   do
   {
      pp1_pow_ui(x, y, nn, p, n, ninv, norm);
      ret = pp1_factor(factor, n, x, nn, norm);
   } while (ret == 1 && factor[0] == 1);

   return ret;
}

int fmpz_factor_pp1(fmpz_t fac, const fmpz_t n_in, ulong B0, ulong c)
{
   long i, j;
   int ret = 0;
   mp_size_t nn = fmpz_size(n_in), r;
   mp_ptr x, y, oldx, oldy, n, ninv, factor;
   ulong pr, oldpr, sqrt, bits0, norm;
   n_primes_t iter;

   if (fmpz_is_even(n_in))
   {
      fmpz_set_ui(fac, 2);
      return 1;
   }

   n_primes_init(iter);

   sqrt = n_sqrt(B0);
   bits0 = FLINT_BIT_COUNT(B0);

   x      = flint_malloc(nn*sizeof(mp_limb_t));
   y      = flint_malloc(nn*sizeof(mp_limb_t));
   oldx   = flint_malloc(nn*sizeof(mp_limb_t));
   oldy   = flint_malloc(nn*sizeof(mp_limb_t));
   n      = flint_malloc(nn*sizeof(mp_limb_t));
   ninv   = flint_malloc(nn*sizeof(mp_limb_t));
   factor = flint_malloc(nn*sizeof(mp_limb_t));
      
   if (nn == 1)
   {
      n[0] = fmpz_get_ui(n_in);
      count_leading_zeros(norm, n[0]);
      n[0] <<= norm;
   } else
   {
      mp_ptr np = COEFF_TO_PTR(*n_in)->_mp_d;
      count_leading_zeros(norm, np[nn - 1]);
      if (norm)
         mpn_lshift(n, np, nn, norm);
      else
         mpn_copyi(n, np, nn);
   }

   flint_mpn_preinvn(ninv, n, nn);
   
   pp1_set_ui(x, nn, norm, c);
   
   /* mul by various prime powers */
   pr = 0;
   oldpr = 0;

   for (i = 0; pr < B0; )
   {
      j = i + 1024;
      oldpr = pr;
      pp1_set(oldx, oldy, x, y, nn);
      for ( ; i < j; i++)
      {
         pr = n_primes_next(iter);
         if (pr < sqrt)
         {
            ulong bits = FLINT_BIT_COUNT(pr);
            ulong exp = bits0 / bits;
            pp1_pow_ui(x, y, nn, n_pow(pr, exp), n, ninv, norm);
         } else
            pp1_pow_ui(x, y, nn, pr, n, ninv, norm);
      }
      
      r = pp1_factor(factor, n, x, nn, norm);
      if (r == 0)
         break;
      if (r != 1 || factor[0] != 1)
      {
         ret = 1;
         goto cleanup;
      }
   }

   if (pr < B0) /* factor = 0 */
   {
      n_primes_jump_after(iter, oldpr);
      pp1_set(x, y, oldx, oldy, nn);
         
      do
      {
         pr = n_primes_next(iter);
         pp1_set(oldx, oldy, x, y, nn);
         if (pr < sqrt)
         {
            ulong bits = FLINT_BIT_COUNT(pr);
            ulong exp = bits0 / bits;
            pp1_pow_ui(x, y, nn, n_pow(pr, exp), n, ninv, norm);
         } else
            pp1_pow_ui(x, y, nn, pr, n, ninv, norm);

         r = pp1_factor(factor, n, x, nn, norm);
         if (r == 0)
            break;
         if (r != 1 || factor[0] != 1)
         {
            ret = (r == 0 ? 0 : 1);
            goto cleanup;
         }
      } while (1);
   } else
      goto cleanup;

   /* factor is still 0 */
   ret = pp1_find_power(factor, oldx, oldy, nn, pr, n, ninv, norm);

cleanup:

   if (ret)
   {
      __mpz_struct * fm = _fmpz_promote(fac);
      mpz_realloc(fm, r);
      mpn_copyi(fm->_mp_d, factor, r);
      fm->_mp_size = r;
   }

   flint_free(x);
   flint_free(y);
   flint_free(oldx);
   flint_free(oldy);
   flint_free(n);
   flint_free(ninv);
   
   n_primes_clear(iter);

   return ret;
}
