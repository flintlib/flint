/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpf_vec.h"

void
_mpf_vec_set(mpf * vec1, const mpf * vec2, slong len)
{
    if (vec1 != vec2)
    {
        slong i;
        for (i = 0; i < len; i++)
            mpf_set(vec1 + i, vec2 + i);
    }
}
