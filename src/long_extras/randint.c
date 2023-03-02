/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <limits.h>

#include "flint.h"
#include "ulong_extras.h"
#include "long_extras.h"

mp_limb_signed_t z_randint(flint_rand_t state, mp_limb_t limit)
{
    mp_limb_signed_t z;

    if ((limit == UWORD(0)) || (limit > WORD_MAX))
        limit = WORD_MAX;

    z = n_randlimb(state) % limit;
    if (n_randint(state, 2))
        z = -z;

    return z;
}
