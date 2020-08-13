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

mp_limb_t
n_mod2_precomp(mp_limb_t a, mp_limb_t n, double npre)
{
    mp_limb_t quot;
    slong rem;

    if (a < n)
        return a;
    if ((mp_limb_signed_t) n < WORD(0))
        return a - n;

    if (n == 1)
    {
        quot = a;
        rem = 0;
    } else
    {
        quot = (mp_limb_t) ((double) a * npre);
        rem  = a - quot * n;
    }
    
    if (rem < (mp_limb_signed_t) (-n))
        quot -= (mp_limb_t) ((double) (-rem) * npre);
    else if (rem >= (slong) n)
        quot += (mp_limb_t) ((double) rem * npre);
    else if (rem < WORD(0))
        return rem + n;
    else
        return rem;
    
    rem = a - quot * n;
    if (rem >= (slong) n)
        return rem - n;
    else if (rem < WORD(0))
        return rem + n;
    else
        return rem;
}
