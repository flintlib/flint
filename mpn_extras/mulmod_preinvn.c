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
