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

#include "d_vec.h"

void
_d_vec_set(double *vec1, const double *vec2, slong len2)
{
    if (vec1 != vec2)
    {
        slong i;
        for (i = 0; i < len2; i++)
            vec1[i] = vec2[i];
    }
}
