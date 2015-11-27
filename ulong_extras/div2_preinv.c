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

    Copyright (C) 2009, 2010, 2015 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

/* 
   Method of Niels Moller and Torbjorn Granlund see paper:
   Improved Division by Invariant Integers: (algorithm 4)
   https://gmplib.org/~tege/division-paper.pdf 
*/

ulong n_div2_preinv(ulong a, ulong n, ulong ninv)
{
   ulong norm, q1, q0, r;

   FLINT_ASSERT(n != 0);

   count_leading_zeros(norm, n);
   n <<= norm;

   {
      const ulong u1 = r_shift(a, FLINT_BITS - norm);
      const ulong u0 = (a << norm);

      umul_ppmm(q1, q0, ninv, u1);
      add_ssaaaa(q1, q0, q1, q0, u1, u0);

      q1++;
      r = u0 - q1*n;

      if (r > q0)
      {
         r += n;
         q1--;
      }

      return (r < n) ? q1 : q1 + 1;
   }
}
