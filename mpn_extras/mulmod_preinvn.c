/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"

/* 
   TODO: speedup mpir's mullow and mulhigh and use instead of mul/mul_n
*/  

void flint_mpn_mulmod_preinvn(mp_ptr r, 
        mp_srcptr a, mp_srcptr b, mp_size_t n, 
        mp_srcptr d, mp_srcptr dinv, ulong norm)
{
   mp_limb_t cy;
   mp_ptr t;
   TMP_INIT;

   TMP_START;
   t = TMP_ALLOC(5*n*sizeof(mp_limb_t));

   if (a == b)
      mpn_sqr(t, a, n);
   else
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

   TMP_END;
}
