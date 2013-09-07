/*=============================================================================

Copyright (C) 2007, 2008 David Harvey (zn_poly)
Copyright (C) 2013 William Hart

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=============================================================================*/

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