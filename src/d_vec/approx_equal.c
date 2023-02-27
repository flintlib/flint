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

int
_d_vec_approx_equal(const double *vec1, const double *vec2, slong len,
                    double eps)
{
    slong i;
    if (vec1 == vec2)
        return 1;

    for (i = 0; i < len; i++)
        if (fabs(vec1[i] - vec2[i]) > eps)
            return 0;

    return 1;
}
