/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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

