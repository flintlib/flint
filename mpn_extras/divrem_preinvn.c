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

    Copyright (C) 2013 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"

void flint_mpn_divrem_n_preinvn(mp_ptr q, mp_ptr r, mp_srcptr a, 
                          mp_srcptr d, mp_size_t n, mp_srcptr dinv)
{
   mp_limb_t cy;
   mp_limb_t ts[150];
   mp_ptr t;

   if (n <= 50)
      t = ts;
   else
      t = flint_malloc(3*n*sizeof(mp_limb_t));

   mpn_mul_n(t + n, a + n, dinv, n);
   mpn_add_n(q, t + 2*n, a + n, n);

   mpn_mul_n(t, q, d, n);
   cy = a[n] - t[n] - mpn_sub_n(r, a, t, n);

   while (cy > 0)
   {
      cy -= mpn_sub_n(r, r, d, n);
      mpn_add_1(q, q, n, 1);
   }
   
   if (mpn_cmp(r, d, n) >= 0)
   {
      mpn_sub_n(r, r, d, n);
      mpn_add_1(q, q, n, 1);
   }

   if (n > 50)
      flint_free(t);
}

mp_limb_t flint_mpn_divrem_preinvn(mp_ptr q, mp_ptr a, mp_size_t m, 
                                      mp_srcptr d, mp_size_t n, mp_srcptr dinv)
{
   mp_limb_t hi = 0;
   mp_ptr t;
   
   if (m <= n)
   {
      if (m == n && mpn_cmp(a, d, n) >= 0)
      {
         mpn_sub_n(a, a, d, n);
         return 1;
      }
      else
      {
         mpn_zero(a + m, n - m);
         return 0;
      }
   } 
   
   if (mpn_cmp(a + m - n, d, n) >= 0)
   {
      mpn_sub_n(a + m - n, a + m - n, d, n);
      hi = 1;
   }
      
   while (m > 2*n)
   {
      flint_mpn_divrem_n_preinvn(q + m - 2*n, a + m - 2*n, a + m - 2*n, d, n, dinv);
      m -= n;
   }

   mpn_zero(a + m, 2*n - m);  /* clear top part in case m != 2n */
      
   if (m != 2*n)
      t = flint_malloc(n*sizeof(mp_limb_t));
   else t = q;

   flint_mpn_divrem_n_preinvn(t, a, a, d, n, dinv);

   if (m != 2*n)
   {
      mpn_copyi(q, t, m - n);
      flint_free(t);
   }

   return hi;
}

