/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_vec.h"

void
_fmpq_vec_set_fmpz_vec(fmpq * res, const fmpz * vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        fmpq_set_fmpz(res + i, vec + i);
}
