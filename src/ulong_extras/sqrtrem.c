/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_sqrtrem(mp_limb_t * r, mp_limb_t a)
{
    mp_limb_t is;

    is = (mp_limb_t) sqrt((double) a);

    is -= (is*is > a);
#if FLINT64
    if (is == UWORD(4294967296)) is--;
#endif
    (*r) = a - is*is;

    return is;
}
