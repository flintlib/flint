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

    Copyright (C) 2015 Vladimir Glazachev
   
******************************************************************************/

#include "aprcl.h"

/* 
    returns smallest integer k satisfies:
        log(n) < (k * (k + 1) * 2^(2 * k)) / (2^(k + 1) - k - 2) + 1
*/        
ulong
_unity_zp_pow_select_k(const fmpz_t n)
{
    ulong bits;
    bits = fmpz_bits(n);

    if (bits <= 8)  return 1;
    if (bits <= 24) return 2;
    if (bits <= 69) return 3;
    if (bits <= 196) return 4;
    if (bits <= 538) return 5;
    if (bits <= 1433) return 6;
    if (bits <= 3714) return 7;
    if (bits <= 9399) return 8;
    if (bits <= 23290) return 9;
    if (bits <= 56651) return 10;
    return 11;
}

