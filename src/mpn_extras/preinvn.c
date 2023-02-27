/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
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
