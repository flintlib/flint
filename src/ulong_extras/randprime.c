/*
    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#define ulong mp_limb_t
#include "flint.h"
#include "ulong_extras.h"


mp_limb_t n_randprime(flint_rand_t state, ulong bits, int proved)
{
    mp_limb_t rand;

    if (bits < 2)
    {
        flint_printf("Exception in n_randprime: attempt to generate prime < 2!\n");
        flint_abort();
    }

    if (bits == FLINT_BITS)
    {
        do { rand = n_randbits(state, bits); }
            while (rand >= UWORD_MAX_PRIME);

        rand = n_nextprime(rand, proved);
    }
    else if (bits == 2)
    {
        rand = 2 + n_randint(state, 2);
    }
    else
    {
        do
        {
            rand = n_randbits(state, bits);
            rand = n_nextprime(rand, proved);
        } while ((rand >> bits) > WORD(0));
    }

    return rand;
}

mp_limb_t n_randtest_prime(flint_rand_t state, int proved)
{
    return n_randprime(state, 2 + n_randint(state, FLINT_BITS - 1), proved);
}
