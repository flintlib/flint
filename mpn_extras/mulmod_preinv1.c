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
#include <gmp.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"

void flint_mpn_mulmod_preinv1(mp_ptr r, 
        mp_srcptr a, mp_srcptr b, mp_size_t n, 
        mp_srcptr d, mp_limb_t dinv, ulong norm)
{
   mp_limb_t q;
   mp_limb_t ts[150];
   mp_ptr t;
   len_t i;

   if (n <= 30)
      t = ts;
   else
      t = flint_malloc(5*n*sizeof(mp_limb_t));

   mpn_mul_n(t, a, b, n);
   if (norm)
      mpn_rshift(t, t, 2*n, norm);

   for (i = 2*n - 1; i >= n; i--)
   {
      flint_mpn_divrem21_preinv(q, t[i], t[i - 1], dinv);
      t[i] -= mpn_submul_1(t + i - n, d, n, q);

      if (mpn_cmp(t + i - n, d, n) >= 0 || t[i] != 0)
      {
         q++;
         t[i] -= mpn_sub_n(t + i - n, t + i - n, d, n);
      }
   }

   mpn_copyi(r, t, n);
}
