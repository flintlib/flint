/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void fmpq_inv(fmpq_t dest, const fmpq_t src)
{
    fmpz tmp;

    if (dest != src)
    {
        fmpq_set(dest, src);
    }

    tmp = dest->num;
    dest->num = dest->den;
    dest->den = tmp;

    if (fmpz_sgn(&dest->den) < 0)
    {
        fmpz_neg(&dest->den, &dest->den);
        fmpz_neg(&dest->num, &dest->num);
    }
}
