/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define MPOLY_INLINES_C

#include <string.h>
#include "mpoly.h"

#undef mpoly_copy_monomials

void mpoly_copy_monomials(ulong * exp1, const ulong * exp2, slong len, slong N)
{
    if (len > 0)
        memcpy(exp1, exp2, N * len * sizeof(ulong));
}
