/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef LONGLONG_MSC_H
#define LONGLONG_MSC_H

#include <stdlib.h>
#include <intrin.h>
#include <immintrin.h>

/* Trailing and leading zeros */
# define flint_clz _CountLeadingZeros64
# define flint_ctz flint_ctz
static inline int flint_ctz(ulong x)
{
    unsigned long int index;
    _BitScanForward64(&index, x);
    return index;
}

/* Byte swap */
# define byte_swap(x) do { (x) = _byteswap_uint64(x); } while (0)

/* Multiplication */
#define umul_ppmm(r1, r0, u, v) \
do \
{ \
    (r0) = (u) * (v); \
    (r1) = __umulh(u, v); \
} while (0)

#define smul_ppmm(r1, r0, u, v) \
do \
{ \
    (r0) = (u) * (v); \
    (r1) = __mulh(u, v); \
} while (0)

#endif
