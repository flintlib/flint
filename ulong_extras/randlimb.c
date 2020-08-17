/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

#if FLINT64

mp_limb_t n_randlimb(flint_rand_t state) 
{   
    state->__randval = (state->__randval*UWORD(13282407956253574709) + UWORD(286824421));
    state->__randval2 = (state->__randval2*UWORD(7557322358563246341) + UWORD(286824421));

    return (state->__randval>>32) + ((state->__randval2>>32) << 32);
}

#else

mp_limb_t n_randlimb(flint_rand_t state) 
{   
    state->__randval = (state->__randval*UWORD(1543932465) +  UWORD(1626832771));
    state->__randval2 = (state->__randval2*UWORD(2495927737) +  UWORD(1626832771));

    return (state->__randval>>16) + ((state->__randval2>>16) << 16);
}

#endif
