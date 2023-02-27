/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"

void
_fmpz_vec_get_mpf_vec(mpf * appv, const fmpz * vec, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        fmpz_get_mpf(appv + i, vec + i);
    }
}
