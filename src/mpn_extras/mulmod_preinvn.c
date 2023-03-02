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

/* 
   TODO: speedup mpir's mullow and mulhigh and use instead of mul/mul_n
*/  

void flint_mpn_mulmod_preinvn(mp_ptr r, 
        mp_srcptr a, mp_srcptr b, mp_size_t n, 
        mp_srcptr d, mp_srcptr dinv, ulong norm)
{
   mp_limb_t cy, p1, p2, b0, b1;
   mp_ptr t;
   TMP_INIT;

   TMP_START;
   t = TMP_ALLOC(5*n*sizeof(mp_limb_t));

   if (n == 2)
   {
      if (norm)
      {
         /* mpn_rshift(b, b, n, norm) */
         b0 = (b[0] >> norm) | (b[1] << (FLINT_BITS - norm));
         b1 = b[1] >> norm;
      } else
      {
         b0 = b[0];
         b1 = b[1];
      }
      
      /* mpn_mul_n(t, a, b, n) */
      umul_ppmm(t[1], t[0], a[0], b0);
      umul_ppmm(t[3], t[2], a[1], b1);
      umul_ppmm(p2, p1, a[0], b1);
      add_sssaaaaaa(t[3], t[2], t[1], t[3], t[2], t[1], 0, p2, p1);
      umul_ppmm(p2, p1, a[1], b0);
      add_sssaaaaaa(t[3], t[2], t[1], t[3], t[2], t[1], 0, p2, p1);
      
      /* mpn_mul_n(t + 3*n, t + n, dinv, n) */
      umul_ppmm(t[7], t[6], t[2], dinv[0]);
      umul_ppmm(t[9], t[8], t[3], dinv[1]);
      umul_ppmm(p2, p1, t[2], dinv[1]);
      add_sssaaaaaa(t[9], t[8], t[7], t[9], t[8], t[7], 0, p2, p1);
      umul_ppmm(p2, p1, t[3], dinv[0]);
      add_sssaaaaaa(t[9], t[8], t[7], t[9], t[8], t[7], 0, p2, p1);

      /* mpn_add_n(t + 4*n, t + 4*n, t + n, n) */
      add_ssaaaa(t[9], t[8], t[9], t[8], t[3], t[2]);

      /* mpn_mul_n(t + 2*n, t + 4*n, d, n) */
      umul_ppmm(t[5], t[4], t[8], d[0]);
      t[6] = t[9]*d[1];
      umul_ppmm(p2, p1, t[8], d[1]);
      add_ssaaaa(t[6], t[5], t[6], t[5], p2, p1);
      umul_ppmm(p2, p1, t[9], d[0]);
      add_ssaaaa(t[6], t[5], t[6], t[5], p2, p1);

      /* cy = t[n] - t[3*n] - mpn_sub_n(r, t, t + 2*n, n) */
      sub_dddmmmsss(cy, r[1], r[0], t[2], t[1], t[0], t[6], t[5], t[4]);

      while (cy > 0)
      {
         /* cy -= mpn_sub_n(r, r, d, n) */
         sub_dddmmmsss(cy, r[1], r[0], cy, r[1], r[0], 0, d[1], d[0]);
      }

      if ((r[1] > d[1]) || (r[1] == d[1] && r[0] >= d[0]))
      {
         /* mpn_sub_n(r, r, d, n) */
         sub_ddmmss(r[1], r[0], r[1], r[0], d[1], d[0]);
      }  
   } else
   {
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
   }

   TMP_END;
}
