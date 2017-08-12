/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2017 Apoorv Mishra

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

/* Current implementation follows */

/*
mp_limb_t n_randint(flint_rand_t state, mp_limb_t limit) 
{
    if (limit == UWORD(0)) return n_randlimb(state);
    else return n_randlimb(state) % limit;
}
*/

/* Current implementation ends */

/* New implementation follows */

mp_limb_t n_randint(flint_rand_t state, mp_limb_t limit) 
{
    if ((limit & (limit - 1)) == 0)
    {
        return n_randlimb(state) & (limit - 1);
    }
    else
    {
        const mp_limb_t rand_max = UWORD_MAX;
        mp_limb_t bucket_size, bucket_num, rand_within_range, temp_rand_max;
        mp_limb_t val1, val2, quotient;
        
        int msb = -1;
        int i;

        temp_rand_max = rand_max;

        for (i = FLINT_BITS - 1; i >= 0; i--)
        {
            if ((UWORD(1)<<i) & limit)
            {
                msb = i;
                break;
            }
        }

        val1 = (UWORD(1)<<(FLINT_BITS - msb - 1));
        val2 = limit*val1;
        quotient = UWORD(0);

        /* First iteration of long division method */
        temp_rand_max -= (val2 - UWORD(1));
        quotient |= val1;
        val1 >>= 1;
        val2 >>= 1;

        quotient |= temp_rand_max/limit;

        bucket_size = quotient;
        bucket_num = bucket_size*limit;
        do
        {
            rand_within_range = n_randlimb(state);
        }
        while (rand_within_range >= bucket_num);

        return rand_within_range/bucket_size;
    }
}

/* New implementation ends */
