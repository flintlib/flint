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

mp_limb_t flint_mpn_preinv1(mp_limb_t d1, mp_limb_t d2)
{
   mp_limb_t q, r[2], p[2], cy;
   
   if (d2 + 1 == 0 && d1 + 1 == 0)
      return 0;

   if (d1 + 1 == 0)
      q = ~d1, r[1] = ~d2;
   else
      udiv_qrnnd(q, r[1], ~d1, ~d2, d1 + 1);

   if (d2 + 1 == 0)
      return q;

   r[0] = 0;

   umul_ppmm(p[1], p[0], q, ~d2);
   cy = mpn_add_n(r, r, p, 2);
 
   p[0] = d2 + 1, p[1] = d1 + (d2 + 1 == 0);
   while (cy || mpn_cmp(r, p, 2) >= 0)
   {
      q++;
      cy -= mpn_sub_n(r, r, p, 2);
   }
   
   return q;
}
