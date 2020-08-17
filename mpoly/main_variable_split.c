/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


void mpoly_main_variable_split_LEX(slong * ind, ulong * pexp, const ulong * Aexp,
             slong l1, slong Alen, const ulong * mults, slong num, slong Abits)
{
    slong i, j = 0, s = 0;
    ulong e, mask = (-UWORD(1)) >> (FLINT_BITS - Abits);

    for (i = 0; i < Alen; i++)
    {
        slong top = Aexp[i] >> (Abits*num);
        while (s < l1 - top)
            ind[s++] = i;
        e = 0;
        for (j = num - 1; j >= 0; j--) {
            e = (e * mults[j]) +  ((Aexp[i] >> (j*Abits)) & mask);
        }
        pexp[i] = e;
    }

    while (s <= l1)
        ind[s++] = Alen;
}

void mpoly_main_variable_split_DEG(slong * ind, ulong * pexp, const ulong * Aexp,
             slong l1, slong Alen, ulong deg, slong num, slong Abits)
{
    slong i, j = 0, s = 0;
    ulong e, mask = (-UWORD(1)) >> (FLINT_BITS - Abits);

    for (i = 0; i < Alen; i++)
    {
        slong top = Aexp[i] >> (Abits*num);
        while (s < l1 - top)
            ind[s++] = i;
        e = 0;
        for (j = num - 1; j >= 1; j--)
            e = (e * deg) + ((Aexp[i] >> (j*Abits)) & mask);
        pexp[i] = e;
    }

    while (s <= l1)
        ind[s++] = Alen;
}

