/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"

/* vec1 = pointwise max of vec1 and vec2 */
void _fmpz_vec_max_inplace(fmpz * vec1, const fmpz * vec2, slong len)
{
    slong i;
    for (i = WORD(0); i < len; i++)
    {
        if (fmpz_cmp(vec1 + i, vec2 + i) < 0)
            fmpz_set(vec1 + i, vec2 + i);
    }
}

/* vec1 = pointwise max of vec2 and vec3 */
void _fmpz_vec_max(fmpz * vec1, const fmpz * vec2, const fmpz * vec3, slong len)
{
    slong i;
    for (i = WORD(0); i < len; i++)
    {
        int cmp = fmpz_cmp(vec2 + i, vec3 + i);
        fmpz_set(vec1 + i, (cmp > 0 ? vec2 : vec3) + i);
    }
}
