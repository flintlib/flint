/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_vec_dot_ptr(fmpz_t c, const fmpz * vec1, fmpz ** const vec2,
                                                       slong offset, slong len)
{
    slong i;

    fmpz_zero(c);
    for (i = 0; i < len; i++)
        fmpz_addmul(c, vec1 + i, vec2[i] + offset);
}
