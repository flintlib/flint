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

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"

void flint_mpn_preinvn(mp_ptr dinv, mp_srcptr d, mp_size_t n)
{
   mp_ptr q, r, d1;
   
   d1 = flint_malloc(n*sizeof(mp_limb_t));
   if (mpn_add_1(d1, d, n, 1)) /* check for d + 1 == 0 */
   {
      mpn_zero(dinv, n);
      flint_free(d1);
      return;
   }

   r = flint_malloc((2*n + 1)*sizeof(mp_limb_t));
   q = flint_malloc((n + 2)*sizeof(mp_limb_t));
 
   mpn_zero(r, 2*n);
   r[2*n] = 1;

   mpn_tdiv_qr(q, r, 0, r, 2*n + 1, d1, n);
   mpn_copyi(dinv, q, n);
   
   flint_free(r);
   flint_free(q);
   flint_free(d1);
}

mp_limb_t flint_mpn_rem2n_preinvn(mp_ptr a, 
                 mp_srcptr b, mp_size_t n, mp_srcptr dinv)
{
   mp_limb_t ret;
   mp_limb_t ts[120];
   mp_ptr t;

   if (n <= 30)
      t = ts;
   else
      t = flint_malloc(4*n*sizeof(mp_limb_t));

   /* ensure { a + m - n, n } < { b, n } */
   if ((ret = (mpn_cmp(a + n, b, n) >= 0)))
      mpn_sub_n(a + n, a + n, b, n);
   
   mpn_mul_n(t, a + n, dinv, n);
   mpn_add_n(t + n, t + n, a + n, n);

   mpn_mul_n(t + 2*n, t + n, b, n);
   mpn_sub_n(a, a, t + 2*n, 2*n);

   while (a[n] > 0 || mpn_cmp(a, b, n) >= 0)
      a[n] -= mpn_sub_n(a, a, b, n);
   
   if (n > 30)
      flint_free(t);

   return ret;
}

void flint_mpn_mulmod_preinvn(mp_ptr r, 
        mp_srcptr a, mp_srcptr b, mp_size_t n, 
        mp_srcptr d, mp_srcptr dinv, ulong norm)
{
   mp_limb_t cy;
   mp_limb_t ts[150];
   mp_ptr t;

   if (n <= 30)
      t = ts;
   else
      t = flint_malloc(5*n*sizeof(mp_limb_t));

   mpn_mul_n(t, a, b, n);
   if (norm)
      mpn_rshift(t, t, 2*n, norm);

   mpn_mul_n(t + 3*n, t + n, dinv, n);
   mpn_add_n(t + 4*n, t + 4*n, t + n, n);

   mpn_mul_n(t + 2*n, t + 4*n, d, n);
   cy = t[n] - t[3*n] - mpn_sub_n(r, t, t + 2*n, n);

   while (cy > 0)
      cy -= mpn_sub_n(r, r, d, n);
   
   if (mpn_cmp(r, d, n) >= 0)
      mpn_sub_n(r, r, d, n);

   if (n > 30)
      flint_free(t);
}

/*void flint_mpn_mulmod_preinvn(mp_ptr r, 
        mp_srcptr a, mp_srcptr b, mp_size_t n, 
        mp_srcptr d, mp_srcptr dinv, ulong norm)
{
   mp_limb_t cy;
   mp_limb_t ts[150];
   mp_ptr t;

   if (n <= 30)
      t = ts;
   else
      t = flint_malloc(5*n*sizeof(mp_limb_t));

   mpn_mul_n(t, a, b, n);
   if (norm)
      mpn_rshift(t, t, 2*n, norm);

   __gmpn_mulhigh_n(t + 3*n, t + n, dinv, n);
   mpn_add_n(t + 4*n, t + 4*n, t + n, n);

   __gmpn_mullow_basecase(t + 2*n, t + 4*n, n, d, n, n + 1);
   cy = t[n] - t[3*n] - mpn_sub_n(r, t, t + 2*n, n);

   while (cy > 0)
      cy -= mpn_sub_n(r, r, d, n);
   
   if (mpn_cmp(r, d, n) >= 0)
      mpn_sub_n(r, r, d, n);

   if (n > 30)
      flint_free(t);
}*/