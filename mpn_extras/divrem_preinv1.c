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

#if !defined(_MSC_VER)
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

mp_limb_t flint_mpn_divrem_preinv1(mp_ptr q, mp_ptr a, mp_size_t m, 
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
      FLINT_ASSERT(n > 0);
      a[i] -= mpn_submul_1(a + i - n, b, n, q[i - n]);

      if (mpn_cmp(a + i - n, b, n) >= 0 || a[i] != 0)
      {
         q[i - n]++;
         a[i] -= mpn_sub_n(a + i - n, a + i - n, b, n);
      }
   }

   return ret;
}

#if !defined(_MSC_VER)
#pragma GCC diagnostic warning "-Wunused-variable"
#endif
