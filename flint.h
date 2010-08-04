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

    Copyright (C) 2009 William Hart

******************************************************************************/

#ifndef FLINT_H
#define FLINT_H

#include <mpir.h>
#include <mpfr.h>
#include "longlong.h"

#define ulong unsigned long

#if __GMP_BITS_PER_MP_LIMB == 64
    #define FLINT_BITS 64
    #define FLINT_D_BITS 53
    #define FLINT64 1
#else 
    #define FLINT_BITS 32
    #define FLINT_D_BITS 31
#endif

#define mp_bitcnt_t unsigned long

/*
  We define this here as there is no mpfr.h
 */
typedef __mpfr_struct mpfr;

#define FLINT_ASSERT(param)

#define FLINT_MAX(x, y) ((x) > (y) ? (x) : (y))
#define FLINT_MIN(x, y) ((x) > (y) ? (y) : (x))
#define FLINT_ABS(x) ((long)(x) < 0 ? (-x) : (x))

#define r_shift(in, shift) \
    ((shift == FLINT_BITS) ? 0L : ((in) >> (shift)))

#define l_shift(in, shift) \
    ((shift == FLINT_BITS) ? 0L : ((in) << (shift)))

static inline 
unsigned int FLINT_BIT_COUNT(mp_limb_t x)
{
   unsigned int zeros = FLINT_BITS;
   if (x) count_leading_zeros(zeros, x);
   return FLINT_BITS - zeros;
}

#if !defined(mpn_zero)
#define mpn_zero(xxx, nnn) \
    do \
    { \
        ulong ixxx; \
        for (ixxx = 0; ixxx < nnn; ixxx++) \
            (xxx)[ixxx] = 0UL; \
    } while (0)
#endif

#define mpn_store(xxx, nnn, yyy) \
   do \
   { \
      ulong ixxx; \
      for (ixxx = 0; ixxx < nnn; ixxx++) \
         (xxx)[ixxx] = yyy; \
   } while (0)

#endif

