/*
    Copyright (C) 2007, 2008 David Harvey (zn_poly)
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
#include "nmod_vec.h"
#include "nmod_poly.h"

/*
   Same as _nmod_poly_KS2_pack(), but requires b <= FLINT_BITS.
*/
void
_nmod_poly_KS2_pack1(mp_ptr res, mp_srcptr op, slong n, slong s,
                ulong b, ulong k, slong r)
{
   /* where to write the next limb */
   mp_ptr dest = res;
   
   /* limb currently being filled */
   mp_limb_t buf;

   /* number of bits used in buf; always in [0, FLINT_BITS) */
   ulong buf_b, buf_b_old;

   /* write leading zero-padding */
   while (k >= FLINT_BITS)
   {
      *dest++ = 0;
      k -= FLINT_BITS;
   }

   buf = 0;

   buf_b = k;
   
   for (; n > 0; n--, op += s)
   {
      /* put low bits of current input into buffer */
      buf += *op << buf_b;
      buf_b_old = buf_b;
      buf_b += b;
      if (buf_b >= FLINT_BITS)
      {
         /* buffer is full; flush it */
         *dest++ = buf;
         buf_b -= FLINT_BITS;
         /* put remaining bits of current input into buffer */
         buf = buf_b_old ? (*op >> (FLINT_BITS - buf_b_old)) : 0;
      }
   }
   
   /* write last limb if it's non-empty */
   if (buf_b)
      *dest++ = buf;

   /* zero-pad up to requested length */
   if (r)
   {
      slong written = dest - res;
      for (; written < r; written++)
         *dest++ = 0;
   }
}

void
_nmod_poly_KS2_pack(mp_ptr res, mp_srcptr op, slong n, slong s,
               ulong b, ulong k, slong r)
{
   /* where to write the next limb */
   mp_ptr dest = res;
   
   /* limb currently being filled */
   mp_limb_t buf;

   /* number of bits used in buf; always in [0, FLINT_BITS) */
   ulong buf_b, buf_b_old;

   if (b <= FLINT_BITS)
   {
      /* use specialised version if b is small enough */
      _nmod_poly_KS2_pack1(res, op, n, s, b, k, r);
      return;
   }
   
   /* write leading zero-padding */
   while (k >= FLINT_BITS)
   {
      *dest++ = 0;
      k -= FLINT_BITS;
   }

   buf = 0;

   buf_b = k;
   
   for (; n > 0; n--, op += s)
   {
      /* put low bits of current input into buffer */
      buf += *op << buf_b;
      buf_b_old = buf_b;
      buf_b += b;
      if (buf_b >= FLINT_BITS)
      {
         /* buffer is full; flush it */
         *dest++ = buf;
         buf_b -= FLINT_BITS;
         /* put remaining bits of current input into buffer */
         buf = buf_b_old ? (*op >> (FLINT_BITS - buf_b_old)) : 0;

         /* write as many extra zeroes as necessary */
         if (buf_b >= FLINT_BITS)
         {
            *dest++ = buf;
            buf = 0;
            buf_b -= FLINT_BITS;
            if (buf_b >= FLINT_BITS)
            {
               *dest++ = 0;
               buf_b -= FLINT_BITS;
           }
         }
      }
   }
   
   /* write last limb if it's non-empty */
   if (buf_b)
      *dest++ = buf;

   /* zero-pad up to requested length */
   if (r)
   {
      slong written = dest - res; 
      for (; written < r; written++)
         *dest++ = 0;
   }
}
