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
   Same as _nmod_poly_KS2_unpack(), but requires b <= FLINT_BITS
   (i.e. writes one word per coefficient)
*/
void
_nmod_poly_KS2_unpack1(mp_ptr res, mp_srcptr op, slong n, ulong b,
                  ulong k)
{
   /* limb we're currently extracting bits from */
   mp_limb_t buf = 0;
   /* number of bits currently in buf; always in [0, FLINT_BITS) */
   ulong buf_b = 0;
   
   /* skip over k leading bits */
   while (k >= FLINT_BITS)
   {
      k -= FLINT_BITS;
      op++;
   }

   if (k)
   {
      buf = *op++;
      buf >>= k;
      buf_b = FLINT_BITS - k;
   }

   if (b == FLINT_BITS)
   {
      /* various special cases */
      if (buf_b)
      {
         for (; n > 0; n--)
         {
            /* we need bits from both sides of a limb boundary */
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + (buf << buf_b);
            buf >>= (FLINT_BITS - buf_b);
         }
      }
      else
      {
         for (; n > 0; n--)
            *res++ = *op++;
      }
   }
   else
   {
      ulong mask = (UWORD(1) << b) - 1;

      for (; n > 0; n--)
      {
         if (b <= buf_b)
         {
            /* buf contains all the bits we need */
            *res++ = buf & mask;
            buf >>= b;
            buf_b -= b;
         }
         else
         {
            /* we need bits from both sides of a limb boundary */
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + ((buf << buf_b) & mask);
            buf >>= (b - buf_b);
            buf_b = FLINT_BITS - (b - buf_b);
         }
      }
   }
}



/*
   Same as _nmod_poly_KS2_unpack(), but requires FLINT_BITS < b <= 2 * FLINT_BITS
   (i.e. writes two words per coefficient)
*/
void
_nmod_poly_KS2_unpack2(mp_ptr res, mp_srcptr op, slong n, ulong b,
                  ulong k)
{
   /* limb we're currently extracting bits from */
   mp_limb_t buf = 0;
   /* number of bits currently in buf; always in [0, FLINT_BITS) */
   ulong buf_b = 0;
   
   /* skip over k leading bits */
   while (k >= FLINT_BITS)
   {
      k -= FLINT_BITS;
      op++;
   }

   if (k)
   {
      buf = *op++;
      buf >>= k;
      buf_b = FLINT_BITS - k;
   }

   if (b == 2 * FLINT_BITS)
   {
      n *= 2;
      
      /* various special cases */
      if (buf_b)
      {
         for (; n > 0; n--)
         {
            /* we need bits from both sides of a limb boundary */
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + (buf << buf_b);
            buf >>= (FLINT_BITS - buf_b);
         }
      }
      else
      {
         for (; n > 0; n--)
            *res++ = *op++;
      }
   }
   else
   {
      ulong mask;
      b -= FLINT_BITS;
      mask = (UWORD(1) << b) - 1;

      for (; n > 0; n--)
      {
         /* shunt one whole limb through first */
         if (buf_b)
         {
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + (buf << buf_b);
            buf >>= (FLINT_BITS - buf_b);
         }
         else
            *res++ = *op++;
      
         /* now handle the fractional limb */
         if (b <= buf_b)
         {
            /* buf contains all the bits we need */
            *res++ = buf & mask;
            buf >>= b;
            buf_b -= b;
         }
         else
         {
            /* we need bits from both sides of a limb boundary */
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + ((buf << buf_b) & mask);
            buf >>= (b - buf_b);
            buf_b = FLINT_BITS - (b - buf_b);
         }
      }
   }
}



/*
   Same as _nmod_poly_KS2_unpack(), but requires 2 * FLINT_BITS < b < 3 * FLINT_BITS
   (i.e. writes three words per coefficient)
*/
void
_nmod_poly_KS2_unpack3(mp_ptr res, mp_srcptr op, slong n, ulong b,
                  ulong k)
{
   /* limb we're currently extracting bits from */
   mp_limb_t buf = 0;
   /* number of bits currently in buf; always in [0, FLINT_BITS) */
   ulong buf_b = 0, mask;

   /* skip over k leading bits */
   while (k >= FLINT_BITS)
   {
      k -= FLINT_BITS;
      op++;
   }

   if (k)
   {
      buf = *op++;
      buf >>= k;
      buf_b = FLINT_BITS - k;
   }

   b -= 2 * FLINT_BITS;
   mask = (UWORD(1) << b) - 1;

   for (; n > 0; n--)
   {
      /* shunt two whole limbs through first */
      if (buf_b)
      {
         ulong temp = buf;
         buf = *op++;
         *res++ = temp + (buf << buf_b);
         buf >>= (FLINT_BITS - buf_b);

         temp = buf;
         buf = *op++;
         *res++ = temp + (buf << buf_b);
         buf >>= (FLINT_BITS - buf_b);
      }
      else
      {
         *res++ = *op++;
         *res++ = *op++;
      }
   
      /* now handle the fractional limb */
      if (b <= buf_b)
      {
         /* buf contains all the bits we need */
         *res++ = buf & mask;
         buf >>= b;
         buf_b -= b;
      }
      else
      {
         /* we need bits from both sides of a limb boundary */
         ulong temp = buf;
         buf = *op++;
         *res++ = temp + ((buf << buf_b) & mask);
         buf >>= (b - buf_b);
         buf_b = FLINT_BITS - (b - buf_b);
      }
   }
}


void
_nmod_poly_KS2_unpack(mp_ptr res, mp_srcptr op, slong n, ulong b,
                 ulong k)
{
   if (b <= FLINT_BITS)
      _nmod_poly_KS2_unpack1 (res, op, n, b, k);
   else if (b <= 2 * FLINT_BITS)
      _nmod_poly_KS2_unpack2 (res, op, n, b, k);
   else    /* b < 3 * FLINT_BITS */
      _nmod_poly_KS2_unpack3 (res, op, n, b, k);
}