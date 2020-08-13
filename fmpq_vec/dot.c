/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpq.h"
#include "fmpq_vec.h"

void
_fmpq_vec_dot(fmpq_t res, const fmpq * vec1, const fmpq * vec2, slong len)
{
    slong i;
    fmpq_zero(res);

    for (i = 0; i < len; i++)
    {
        fmpq_addmul(res, vec1 + i, vec2 + i);
    }
}
