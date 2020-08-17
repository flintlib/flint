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

mp_limb_t n_mod_precomp(mp_limb_t a, mp_limb_t n, double npre)
{
    mp_limb_t quot, rem;

    quot = (mp_limb_t) ((double) a * npre);
    rem  = a - quot*n;
    if ((slong) rem < 0) /* unlikely */
       rem += n;
    return rem - (n & (((mp_limb_signed_t) (n - rem - 1)) >> (FLINT_BITS-1)));
}
