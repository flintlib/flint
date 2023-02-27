/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2017 Apoorv Mishra

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

ulong n_randint(flint_rand_t state, ulong limit) 
{
    if (limit == UWORD(0)) return n_randlimb(state);
    else return n_randlimb(state) % limit;
}

mp_limb_t n_urandint(flint_rand_t state, mp_limb_t limit) 
{
    if ((limit & (limit - 1)) == 0)
    {
        return n_randlimb(state) & (limit - 1);
    }
    else
    {
        const mp_limb_t rand_max = UWORD_MAX;
        mp_limb_t bucket_size, num_of_buckets, rand_within_range;

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
