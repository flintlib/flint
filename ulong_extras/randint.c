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
    mp_limb_t rand_max = UWORD_MAX;

    if (limit == UWORD(0)) 
    {
        return n_randlimb(state);
    } 
    else
    {
        mp_limb_t x, y, r;

        x = rand_max/limit;
        y = x*limit;

        do 
        {
            r = n_randlimb(state);
        }
        while (r >= y);

        return r/x;
    }
}

/* New implementation ends */
