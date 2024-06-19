/*
    Copyright (C) 2012 William Hart
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
TODO:
 * fixed-length code for more small n
 * note: 1x1 cannot use nmod_mul because the inverses are defined
         differently
 * use mullow
*/

void flint_mpn_mulmod_preinvn(mp_ptr r,
        mp_srcptr a, mp_srcptr b, mp_size_t n,
        mp_srcptr d, mp_srcptr dinv, ulong norm)
{
   mp_limb_t cy, b0, b1;

   if (n == 2)
   {
      mp_limb_t t[10];

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
      FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], a[1], a[0], b1, b0);

      /* mpn_mul_n(t + 3*n, t + n, dinv, n) */
      FLINT_MPN_MUL_2X2(t[9], t[8], t[7], t[6], t[3], t[2], dinv[1], dinv[0]);

      /* mpn_add_n(t + 4*n, t + 4*n, t + n, n) */
      add_ssaaaa(t[9], t[8], t[9], t[8], t[3], t[2]);

      /* mpn_mul_n(t + 2*n, t + 4*n, d, n) */
      FLINT_MPN_MUL_3P2X2(t[6], t[5], t[4], t[9], t[8], d[1], d[0]);

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
     mp_ptr t;
     TMP_INIT;

     TMP_START;
     t = TMP_ALLOC(5*n*sizeof(mp_limb_t));

      if (a == b)
         flint_mpn_sqr(t, a, n);
      else
         flint_mpn_mul_n(t, a, b, n);

      if (norm)
         mpn_rshift(t, t, 2*n, norm);

      flint_mpn_mul_or_mulhigh_n(t + 3*n, t + n, dinv, n);
      mpn_add_n(t + 4*n, t + 4*n, t + n, n);

      /* note: we rely on the fact that mul_or_mullow_n actually
               writes at least n + 1 limbs */
      flint_mpn_mul_or_mullow_n(t + 2*n, t + 4*n, d, n);
      cy = t[n] - t[3*n] - mpn_sub_n(r, t, t + 2*n, n);

      while (cy > 0)
         cy -= mpn_sub_n(r, r, d, n);

      if (mpn_cmp(r, d, n) >= 0)
         mpn_sub_n(r, r, d, n);

      FLINT_ASSERT(mpn_cmp(r, d, n) < 0);

     TMP_END;
   }
}

void flint_mpn_mulmod_preinvn_2(mp_ptr r,
        mp_srcptr a, mp_srcptr b,
        mp_srcptr d, mp_srcptr dinv, ulong norm)
{
    mp_limb_t cy, b0, b1, r0, r1;
    mp_limb_t t[10];

    if (norm)
    {
        /* mpn_lshift(b, b, n, norm) */
        b0 = (b[0] << norm);
        b1 = (b[1] << norm) | (b[0] >> (FLINT_BITS - norm));
    }
    else
    {
        b0 = b[0];
        b1 = b[1];
    }

    /* mpn_mul_n(t, a, b, n) */
    FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], a[1], a[0], b1, b0);

    /* mpn_mul_n(t + 3*n, t + n, dinv, n) */
    FLINT_MPN_MUL_2X2(t[9], t[8], t[7], t[6], t[3], t[2], dinv[1], dinv[0]);

    /* mpn_add_n(t + 4*n, t + 4*n, t + n, n) */
    add_ssaaaa(t[9], t[8], t[9], t[8], t[3], t[2]);

    /* mpn_mul_n(t + 2*n, t + 4*n, d, n) */
    FLINT_MPN_MUL_3P2X2(t[6], t[5], t[4], t[9], t[8], d[1], d[0]);

    /* cy = t[n] - t[3*n] - mpn_sub_n(r, t, t + 2*n, n) */
    sub_dddmmmsss(cy, r1, r0, t[2], t[1], t[0], t[6], t[5], t[4]);

    while (cy > 0)
    {
        /* cy -= mpn_sub_n(r, r, d, n) */
        sub_dddmmmsss(cy, r1, r0, cy, r1, r0, 0, d[1], d[0]);
    }

    if ((r1 > d[1]) || (r1 == d[1] && r0 >= d[0]))
    {
        /* mpn_sub_n(r, r, d, n) */
        sub_ddmmss(r1, r0, r1, r0, d[1], d[0]);
    }

    if (norm)
    {
        r[0] = (r0 >> norm) | (r1 << (FLINT_BITS - norm));
        r[1] = (r1 >> norm);
    }
    else
    {
        r[0] = r0;
        r[1] = r1;
    }
}
