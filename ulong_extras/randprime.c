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

    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#undef ulong /* prevent clash with standard library */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long
#include "flint.h"
#include "ulong_extras.h"


mp_limb_t n_randprime(flint_rand_t state, unsigned long bits, int proved)
{
    mp_limb_t rand;

    if (bits < 2)
    {
        printf("Exception in n_randprime: attempt to generate prime < 2!\n");
        abort();
    }

    if (bits == FLINT_BITS)
    {
        do { rand = n_randbits(state, bits); }
            while (rand >= ULONG_MAX_PRIME);

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
        } while ((rand >> bits) > 0L);
    }

    return rand;
}

mp_limb_t n_randtest_prime(flint_rand_t state, int proved)
{
    return n_randprime(state, 2 + n_randint(state, FLINT_BITS - 1), proved);
}
