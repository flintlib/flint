/*
    Copyright (C) 2008, 2009, 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

int
_fmpz_vec_equal(const fmpz * vec1, const fmpz * vec2, slong len)
{
    slong i;
    if (vec1 == vec2)
        return 1;

    for (i = 0; i < len; i++)
        if (!fmpz_equal(vec1 + i, vec2 + i))
            return 0;

    return 1;
}
