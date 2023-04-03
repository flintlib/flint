/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
fmpq_height(fmpz_t height, const fmpq_t x)
{
    if (fmpz_cmpabs(fmpq_numref(x), fmpq_denref(x)) < 0)
        fmpz_abs(height, fmpq_denref(x));
    else
        fmpz_abs(height, fmpq_numref(x));
}
