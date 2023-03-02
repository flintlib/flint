/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "mpoly.h"

void mpoly_pack_monomials_tight(ulong * exp1, const ulong * exp2,
                     slong len, const slong * mults, slong nfields, slong bits)
{
    slong i, j, shift;
    ulong e, mask = (-UWORD(1)) >> (FLINT_BITS - bits);

    for (i = 0; i < len; i++)
    {
        shift = (nfields - 1)*bits;
        e = (exp2[i] >> shift) & mask;
        for (j = nfields - 2; j >= 0; j--)
        {
            shift -= bits;
            e *= mults[j];
            e += (exp2[i] >> shift) & mask;
        }
        exp1[i] = e;
   }
}


