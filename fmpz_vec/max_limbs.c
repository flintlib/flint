/*
    Copyright (C) 2010 William Hart

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

mp_size_t
_fmpz_vec_max_limbs(const fmpz * vec, slong len)
{
    slong i;
    mp_size_t limbs, max_limbs = 0;

    for (i = 0; i < len; i++)
    {
        limbs = fmpz_size(vec + i);
        if (limbs > max_limbs)
            max_limbs = limbs;
    }

    return max_limbs;
}
