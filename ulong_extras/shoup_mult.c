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

    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

/* Computes the W' = [w * b / p] (b = mp_limb_t power) */
mp_limb_t
n_mulmod_precomp_shoup(mp_limb_t w, mp_limb_t p)
{
   mp_limb_t q, r;
   udiv_qrnnd(q, r, w, UWORD(0), p);
   return q;
}

/* Computes the w * t mod p.
   w, t < p; p < (mp_limb_t power) / 2. 
*/
mp_limb_t
n_mulmod_shoup(mp_limb_t w, mp_limb_t t, mp_limb_t w_precomp, mp_limb_t p)
{
   mp_limb_t q, r, p_hi, p_lo;

   umul_ppmm(p_hi, p_lo, w_precomp, t);
   q = p_hi;

   
   r = w * t;
   r -= q * p;

   if (r >= p)
   {
      r -= p;
   }

    return r;
}
