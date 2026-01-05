/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

void radix_init(radix_t radix, ulong b, unsigned int exp)
{
    ulong B;
    ulong hi;
    int i;

    if (b < 2 || exp >= FLINT_BITS)
        flint_throw(FLINT_ERROR, "radix_init: require b >= 2 and exp < FLINT_BITS");

    B = b;
    for (i = 2; i <= exp; i++)
    {
        umul_ppmm(hi, B, B, b);
        if (hi != 0)
            flint_throw(FLINT_ERROR, "radix_init: require b^e < 2^FLINT_BITS");
    }

    nmod_init(&radix->b, b);
    radix->exp = exp;
    nmod_init(&radix->B, B);
}

void radix_clear(radix_t radix)
{
}
