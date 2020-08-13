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

mp_limb_t n_mulmod_precomp(mp_limb_t a, mp_limb_t b, mp_limb_t n, double npre)
{
    mp_limb_t quot;
    mp_limb_signed_t rem;

    quot = (mp_limb_t) ((double) a * (double) b * npre);
    rem  = a * b - quot * n;
    if (rem < 0) 
    {
        rem += n;
        if (rem < 0) return rem + n;
    }
    else if (rem >= n) return rem - n;
    return rem;
}
