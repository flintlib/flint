/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"

#if !defined(_MSC_VER)
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

void flint_mpn_mulmod_preinv1(mp_ptr r,
        mp_srcptr a, mp_srcptr b, mp_size_t n,
        mp_srcptr d, mp_limb_t dinv, ulong norm)
{
   mp_limb_t q;
   mp_limb_t ts[150];
   mp_ptr t;
   slong i;

   FLINT_ASSERT(n > 0);

   if (n <= 30)
      t = ts;
   else
      t = flint_malloc(5*n*sizeof(mp_limb_t));

   if (a == b)
      flint_mpn_sqr(t, a, n);
   else
      flint_mpn_mul_n(t, a, b, n);

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

   if (n > 30)
       flint_free(t);
}

#if !defined(_MSC_VER)
#pragma GCC diagnostic warning "-Wunused-variable"
#endif
