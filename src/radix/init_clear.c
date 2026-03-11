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
    ulong hi, lo;
    int i;

    if (b < 2 || exp >= FLINT_BITS)
        flint_throw(FLINT_ERROR, "radix_init: require b >= 2 and exp < FLINT_BITS");

    B = b;

    if (exp == 0)
    {
        for (i = 1; ; i++)
        {
            umul_ppmm(hi, lo, B, b);
            if (hi != 0)
            {
                exp = i;
                break;
            }
            else
            {
                B = lo;
            }
        }
    }
    else
    {
        for (i = 2; i <= exp; i++)
        {
            umul_ppmm(hi, B, B, b);
            if (hi != 0)
                flint_throw(FLINT_ERROR, "radix_init: require b^e < 2^FLINT_BITS");
        }
    }

    nmod_init(&radix->b, b);
    radix->exp = exp;
    nmod_init(&radix->B, B);

    radix->bpow = flint_malloc(sizeof(ulong) * (exp + 1));
    radix->bpow_div = flint_malloc(sizeof(n_div_precomp_struct) * (exp + 1));

    /* todo: combine with earlier power computations */
    radix->bpow[0] = 1;
    radix->bpow[1] = b;
    for (i = 2; i <= exp; i++)
        radix->bpow[i] = radix->bpow[i - 1] * b;

    for (i = 0; i <= exp; i++)
        n_div_precomp_init(radix->bpow_div + i, radix->bpow[i]);


    {
        int prevbc = -1;
        for (i = 0; i <= exp; i++)
        {
            unsigned int bc = FLINT_BIT_COUNT(radix->bpow[i]);
            int jbc;

            for (jbc = prevbc + 1; jbc <= bc; jbc++)
                radix->bits_to_digit_size[jbc - 1] = i;

            prevbc = bc;
        }
    }
}

void radix_clear(radix_t radix)
{
    flint_free(radix->bpow);
    flint_free(radix->bpow_div);
}

