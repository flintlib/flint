/* 
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fft.h"
#include "longlong.h"

/* WARNING: relies on GCC's handling of >> as arithmetic shift right */

void mpn_div_2expmod_2expp1(ulong * t, ulong * i1, mp_size_t limbs, flint_bitcnt_t d)
{
   ulong lo;
   ulong * ptr;
   mp_limb_signed_t hi;
   
   if (d == 0)
   {   
      if (t != i1)
         FLINT_MPN_COPYI(t, i1, limbs + 1);
   } else
   {
      hi = i1[limbs];
      lo = mpn_rshift(t, i1, limbs + 1, d);
      t[limbs] = (hi >> d);
      ptr = t + limbs - 1;
      sub_ddmmss(ptr[1], ptr[0], ptr[1], ptr[0], (ulong) 0, lo);
   }
}
