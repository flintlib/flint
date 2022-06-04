/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec-impl.h"


void
_fmpz_vec_scalar_tdiv_q_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_tdiv_q_si(vec1 + i, vec2 + i, c);
}
