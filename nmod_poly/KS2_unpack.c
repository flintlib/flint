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
