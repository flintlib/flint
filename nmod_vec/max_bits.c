/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"

mp_bitcnt_t _nmod_vec_max_bits(mp_srcptr vec, slong len)
{
    mp_bitcnt_t bits = 0;
    mp_limb_t mask   = ~(mp_limb_t) 0;
    slong i;

    for (i = 0; i < len; i++)
    {
        if (vec[i] & mask)
        {
            bits = FLINT_BIT_COUNT(vec[i]);
            if (bits == FLINT_BITS) break;
            else mask = ~(mp_limb_t) 0 - ((UWORD(1) << bits) - UWORD(1));
        }
    }

    return bits;
}
