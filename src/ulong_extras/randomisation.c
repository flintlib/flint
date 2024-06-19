/*
    Copyright (C) 2007-2009 William Hart
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2017 Apoorv Mishra

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"

#define l_shift(in, shift) \
    ((shift == FLINT_BITS) ? WORD(0) : ((in) << (shift)))

ulong n_randbits(flint_rand_t state, unsigned int bits)
{
   if (bits == 0)
       return UWORD(0);
   else
       return (UWORD(1) << (bits - 1)) | n_randint(state, l_shift(UWORD(1), bits));
}

ulong n_urandint(flint_rand_t state, ulong limit)
{
    if ((limit & (limit - 1)) == 0)
    {
        return n_randlimb(state) & (limit - 1);
    }
    else
    {
        const ulong rand_max = UWORD_MAX;
        ulong bucket_size, num_of_buckets, rand_within_range;

        bucket_size = 1 + (rand_max - limit + 1)/limit;
        num_of_buckets = bucket_size*limit;
        do
        {
            rand_within_range = n_randlimb(state);
        }
        while (rand_within_range >= num_of_buckets);

        return rand_within_range/bucket_size;
    }
}

#if FLINT64
ulong n_randlimb(flint_rand_t state)
{
    state->__randval = (state->__randval*UWORD(13282407956253574709) + UWORD(286824421));
    state->__randval2 = (state->__randval2*UWORD(7557322358563246341) + UWORD(286824421));

    return (state->__randval>>32) + ((state->__randval2>>32) << 32);
}
#else
ulong n_randlimb(flint_rand_t state)
{
    state->__randval = (state->__randval*UWORD(1543932465) +  UWORD(1626832771));
    state->__randval2 = (state->__randval2*UWORD(2495927737) +  UWORD(1626832771));

    return (state->__randval>>16) + ((state->__randval2>>16) << 16);
}
#endif

ulong n_randtest_bits(flint_rand_t state, int bits)
{
    ulong m;
    ulong n;

    m = n_randlimb(state);

    if (m & UWORD(7))
    {
        n = n_randbits(state, bits);
    }
    else
    {
        m >>= 3;

        switch (m & UWORD(7))
        {
            case 0:  n = 0;         break;
            case 1:  n = 1;         break;
            case 2:  n = COEFF_MAX; break;
            case 3:  n = WORD_MAX;  break;
            case 4:  n = UWORD_MAX; break;
            case 5:  n =  (UWORD(1)<<n_randint(state, FLINT_BITS))
                        - (UWORD(1)<<n_randint(state, FLINT_BITS));
                                    break;
            case 6:  n =  (UWORD(1)<<n_randint(state, FLINT_BITS));
                                    break;
            case 7:  n = -(UWORD(1)<<n_randint(state, FLINT_BITS));
                                    break;
            default: n = 0;
        }

        if (bits < FLINT_BITS)
           n &= ((UWORD(1)<<bits) - UWORD(1)); /* mask it off */

        if (bits) /* set most significant bit */
           n |= (UWORD(1)<<(bits - 1));
        else
           n = 0;
    }

    return n;
}

ulong n_randtest(flint_rand_t state)
{
    return n_randtest_bits(state, n_randint(state, FLINT_BITS + 1));
}

ulong n_randtest_not_zero(flint_rand_t state)
{
    ulong n;

    while ((n = n_randtest(state)) == 0) ;
    return n;
}

ulong n_randprime(flint_rand_t state, ulong bits, int proved)
{
    ulong rand;

    if (bits < 2)
    {
        flint_throw(FLINT_ERROR, "Exception in n_randprime: attempt to generate prime < 2!\n");
    }

    if (bits == FLINT_BITS)
    {
        do
            rand = n_randbits(state, bits);
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

ulong n_randtest_prime(flint_rand_t state, int proved)
{
    return n_randprime(state, 2 + n_randint(state, FLINT_BITS - 1), proved);
}
