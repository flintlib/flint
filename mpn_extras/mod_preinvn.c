/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# else
#  if HAVE_ALLOCA_H
#   include <alloca.h>
#  else
#   if _MSC_VER
#    include <malloc.h>
#    define alloca _alloca
#   else
#    ifdef __DECC
#     define alloca(x) __ALLOCA(x)
#    else
#     ifdef BSD
#      include <stdlib.h>
#     else
#      error Could not find alloca
#     endif
#    endif
#   endif
#  endif
# endif
#endif

#include "mpn_extras.h"
#include "flint-impl.h"

/* 
   TODO: speedup mpir's mullow and mulhigh and use instead of mul/mul_n
*/  

void flint_mpn_mod_preinvn(ulong_ptr rp, ulong_srcptr ap, mp_size_t m, 
                          ulong_srcptr d, mp_size_t n, ulong_srcptr dinv)
{
   ulong cy;
   ulong_ptr t, r, a;
   mp_size_t size;
   TMP_INIT;

   a = (ulong_ptr) ap + m - 2*n;
   r = rp + m - 2*n;

   /* check if top n limbs of a exceed d */
   if (mpn_cmp(a + n, d, n) >= 0)
      mpn_sub_n(r + n, a + n, d, n);
   else if (r != a)
      mpn_copyi(r + n, a + n, n);

   TMP_START;
   t = TMP_ALLOC(3*n*sizeof(ulong));

   /* 2n by n division */
   while (m >= 2*n)
   {
      mpn_mul_n(t, dinv, r + n, n);
      cy = mpn_add_n(t + 2*n, t + n, r + n, n);

      mpn_mul_n(t, d, t + 2*n, n);
      cy = r[n] - t[n] - mpn_sub_n(r, a, t, n);

      while (cy > 0)
         cy -= mpn_sub_n(r, r, d, n);

      if (mpn_cmp(r, d, n) >= 0)
         mpn_sub_n(r, r, d, n);

      m -= n;
      r -= n;
      a -= n;
   }

   size = m - n;

   /* m by n division with 2n > m > n */
   if (size)
   {
      if (rp != ap)
         mpn_copyi(rp, ap, size);
      
      mpn_mul(t, dinv, n, rp + n, size);
      cy = mpn_add_n(t + 2*n, t + n, rp + n, size);

      mpn_mul(t, d, n, t + 2*n, size);
      if (cy)
         mpn_add_n(t + size, t + size, d, n + 1 - size);
      
      cy = rp[n] - t[n] - mpn_sub_n(rp, rp, t, n);

      while (cy > 0)
         cy -= mpn_sub_n(rp, rp, d, n);

      if (mpn_cmp(rp, d, n) >= 0)
         mpn_sub_n(rp, rp, d, n);
   }

   TMP_END;
}
