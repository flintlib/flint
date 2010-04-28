/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2009 William Hart

*****************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

const ulong revtab0[1] = { 0 };
const ulong revtab1[2] = { 0, 1 };
const ulong revtab2[4] = { 0, 2, 1, 3 };
const ulong revtab3[8] = { 0, 4, 2, 6, 1, 5, 3, 7 };
const ulong revtab4[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

const ulong * revtab[5] = { revtab0, revtab1, revtab2, revtab3, revtab4 };

// computes the reverse binary of a binary number of the given number of bits
ulong n_revbin(ulong in, ulong bits)
{
   ulong out = 0, i;
   
   if (bits <= 4)
      return revtab[bits][in];

   for (i = 0; i < bits; i++)
   {   
      out <<= 1;
      out += (in & 1);
      in >>= 1;
   }

   return out;
}
