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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"

mp_bitcnt_t _nmod_vec_max_bits(mp_srcptr vec, len_t len)
{
   mp_bitcnt_t bits = 0;
   mp_limb_t mask   = ~(mp_limb_t) 0;
   len_t i;
   
   for (i = 0; i < len; i++)
   {
      if (vec[i] & mask)
      {
         bits = FLINT_BIT_COUNT(vec[i]);
         if (bits == FLINT_BITS) break;
         else mask = ~(mp_limb_t) 0 - ((1UL << bits) - 1UL);
      }
   }
   
   return bits;
}
