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

#include "flint.h"

int
_mpf_vec_approx_equal(mpf_srcptr vec1, mpf_srcptr vec2, slong len,
                      flint_bitcnt_t bits)
{
    slong i;
    if (vec1 == vec2)
        return 1;

    for (i = 0; i < len; i++)
        if (mpf_eq(vec1 + i, vec2 + i, bits) == 0)
            return 0;

    return 1;
}
