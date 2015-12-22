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

    Copyright (C) 2009, 2015 William Hart
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

#if (defined (__amd64__) || defined (__i386__) || defined (__i486__))	 

ulong n_gcd(ulong x, ulong y)
{
   register ulong s0, s1, f;

   if (x == 0) return y;
   if (y == 0) return x;

   count_trailing_zeros(s0, x);
   count_trailing_zeros(s1, y);

   f = FLINT_MIN(s0, s1);

   x >>= s0;
   y >>= s1;

   while (x != y)
   {
      if (x < y)
      {
         y -= x;
         count_trailing_zeros(s1, y);
         y >>= s1;
      } else
      {
         x -= y;
         count_trailing_zeros(s0, x);
         x >>= s0;
      }
   }

   return x << f;
}

#else

ulong n_gcd(ulong x, ulong y)
{
   ulong d, v, quot, rem;

   if (x >= y)
      v = y;
   else
   {
      v = x;
      x = y;
   }

   /* x and y both have top bit set */
   if ((slong) (x & v) < WORD(0))  
   {
      d = x - v;
      x = v;
      v = d;
   }

   /* second value has second msb set */
   while ((slong) (v << 1) < WORD(0))  
   {
      d = x - v;
      x = v;
      if (d < v)                v = d;            /* quot = 1 */
      else if (d < (v << 1))    v = d - x;        /* quot = 2 */
      else                      v = d - (x << 1); /* quot = 3 */
   }

   while (v)
   {
      /* overflow not possible due to top 2 bits of v not being set */
      if (x < (v << 2)) /* avoid divisions when quotient < 4 */
      {
         d = x - v;
         x = v;
         if (d < v)             v = d;            /* quot = 1 */
         else if (d < (v << 1)) v = d - x;        /* quot = 2 */
         else                   v = d - (x << 1); /* quot = 3 */
      } else
      {
         quot = x / v;
         rem  = x - v * quot;
         x    = v;
         v    = rem;
      }
   }

   return x;
}

#endif
