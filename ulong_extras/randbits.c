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

mp_limb_t n_randbits(flint_rand_t state, unsigned int bits)
{
   if (bits == 0) return UWORD(0);
   else return (UWORD(1) << (bits - 1)) | n_randint(state, l_shift(UWORD(1), bits));
}
