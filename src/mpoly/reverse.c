/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


void mpoly_reverse(ulong * Aexp, const ulong * Bexp, slong len, slong N)
{
    slong i;

    if (Aexp == Bexp)
    {
        for (i = 0; i < len/2; i++)
        {
            mpoly_monomial_swap(Aexp + N*i, Aexp + N*(len - i - 1), N);
        }
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            mpoly_monomial_set(Aexp + N*i, Bexp + N*(len - i - 1), N);
        }
    }
}
