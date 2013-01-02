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

mp_limb_t flint_mpn_preinv1(mp_limb_t d1, mp_limb_t d2)
{
   mp_limb_t q, r[2], p[2], cy;
   
   if (d2 + 1 == 0 && d1 + 1 == 0)
      return 0;

   if (d1 + 1 == 0)
      q = ~d1, r[1] = ~d2;
   else
      udiv_qrnnd(q, r[1], ~d1, ~d2, d1 + 1);

   r[0] = 0;

   if (d2 + 1 == 0)
      add_ssaaaa(cy, r[1], 0, r[1], 0, q);   
   else
   {
      umul_ppmm(p[1], p[0], q, ~d2 - 1);
      cy = mpn_add_n(r, r, p, 2);
   }
 
   p[0] = d2 + 1, p[1] = d1 + (d2 + 1 == 0);
   if (cy || mpn_cmp(r, p, 2) >= 0)
      q++;
   
   return q;
}

mp_limb_t flint_mpn_divrem_basecase_preinv1(mp_ptr q, mp_ptr a, mp_size_t m, 
                                 mp_srcptr b, mp_size_t n, mp_limb_t dinv)
{
   mp_limb_t ret;
   mp_size_t i;

   /* ensure { a + i, n } < { b, n } */
   if ((ret = (mpn_cmp(a + m - n, b, n) >= 0)))
      mpn_sub_n(a + m - n, a + m - n, b, n);
   
   for (i = m - 1; i >= n; i--)
   {
      flint_mpn_divrem21_preinv(q[i - n], a[i], a[i - 1], dinv);
      a[i] -= mpn_submul_1(a + i - n, b, n, q[i - n]);

      if (mpn_cmp(a + i - n, b, n) >= 0 || a[i] != 0)
      {
         q[i - n]++;
         a[i] -= mpn_sub_n(a + i - n, a + i - n, b, n);
      }
   }

   return ret;
}