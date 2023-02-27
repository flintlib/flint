/*
    Copyright (C) 2013 William Hart

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

/* 
   TODO: speedup mpir's mullow and mulhigh and use instead of mul/mul_n
*/  

mp_limb_t flint_mpn_divrem_preinvn(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t m,
                          mp_srcptr d, mp_size_t n, mp_srcptr dinv)
{
   mp_limb_t cy, hi = 0;
   mp_ptr t, q, r, a;
   mp_size_t size;
   TMP_INIT;

   a = (mp_ptr) ap + m - 2*n;
   r = rp + m - 2*n;

   /* check if top n limbs of a exceed d */
   if (mpn_cmp(a + n, d, n) >= 0)
   {
      mpn_sub_n(r + n, a + n, d, n);
      hi = 1;
   } else if (r != a)
      mpn_copyi(r + n, a + n, n);

   q = qp + m - 2*n;

   TMP_START;
   t = TMP_ALLOC(2*n*sizeof(mp_limb_t));

   /* 2n by n division */
   while (m >= 2*n)
   {
      mpn_mul_n(t, dinv, r + n, n);
      cy = mpn_add_n(q, t + n, r + n, n);

      mpn_mul_n(t, d, q, n);
      cy = r[n] - t[n] - mpn_sub_n(r, a, t, n);

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

      m -= n;
      r -= n;
      a -= n;
      q -= n;
   }

   size = m - n;

   /* m by n division with 2n > m > n */
   if (size)
   {
      if (rp != ap)
         mpn_copyi(rp, ap, size);
      
      mpn_mul(t, dinv, n, rp + n, size);
      cy = mpn_add_n(qp, t + n, rp + n, size);

      mpn_mul(t, d, n, qp, size);
      if (cy)
         mpn_add_n(t + size, t + size, d, n + 1 - size);
      
      cy = rp[n] - t[n] - mpn_sub_n(rp, rp, t, n);

      while (cy > 0)
      {
         cy -= mpn_sub_n(rp, rp, d, n);
         mpn_add_1(qp, qp, size, 1);
      }

      if (mpn_cmp(rp, d, n) >= 0)
      {
         mpn_sub_n(rp, rp, d, n);
         mpn_add_1(qp, qp, size, 1);
      }
   }

   TMP_END;

   return hi;
}
