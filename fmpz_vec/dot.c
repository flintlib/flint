/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"

void
_fmpz_vec_dot(fmpz_t res, const fmpz * vec1, const fmpz * vec2, slong len2)
{
    slong i;

    fmpz_zero(res);
    for (i = 0; i < len2; i++)
        fmpz_addmul(res, vec1 + i, vec2 + i);
}
