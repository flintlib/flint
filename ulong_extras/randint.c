/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_randint(flint_rand_t state, mp_limb_t limit) 
{
   if (limit == UWORD(0)) return n_randlimb(state);
   else return n_randlimb(state) % limit;
}
