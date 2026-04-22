/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "ulong_extras.h"

static const uint8_t n_sizeinbase10_tab1[64] = {
    1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8,
    8, 8, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 13,
    14, 14, 14, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19,
    19, 19, 20,
};

static const ulong n_sizeinbase10_tab2[FLINT_BITS] = {
    1, 1, 1, 10, 10, 10, 100, 100, 100, 1000, 1000, 1000, 1000, 10000, 10000,
    10000, 100000, 100000, 100000, 1000000, 1000000, 1000000, 1000000, 10000000,
    10000000, 10000000, 100000000, 100000000, 100000000, 1000000000,
    1000000000, 1000000000,
#if FLINT_BITS == 64
    1000000000, UWORD(10000000000), UWORD(10000000000), UWORD(10000000000),
    UWORD(100000000000), UWORD(100000000000), UWORD(100000000000),
    UWORD(1000000000000), UWORD(1000000000000), UWORD(1000000000000),
    UWORD(1000000000000), UWORD(10000000000000), UWORD(10000000000000),
    UWORD(10000000000000), UWORD(100000000000000), UWORD(100000000000000),
    UWORD(100000000000000), UWORD(1000000000000000), UWORD(1000000000000000),
    UWORD(1000000000000000), UWORD(1000000000000000), UWORD(10000000000000000),
    UWORD(10000000000000000), UWORD(10000000000000000),
    UWORD(100000000000000000), UWORD(100000000000000000),
    UWORD(100000000000000000),  UWORD(1000000000000000000),
    UWORD(1000000000000000000), UWORD(1000000000000000000),
    UWORD(1000000000000000000), UWORD(10000000000000000000),
#endif
};

slong
n_nonzero_sizeinbase10(ulong n)
{
    int b;
    FLINT_ASSERT(n != 0);
    b = FLINT_BIT_COUNT(n) - 1;
    return n_sizeinbase10_tab1[b] - (n < n_sizeinbase10_tab2[b]);
}
