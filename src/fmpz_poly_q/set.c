/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_q.h"

void fmpz_poly_q_set(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    if (rop != op)
    {
        fmpz_poly_set(rop->num, op->num);
        fmpz_poly_set(rop->den, op->den);
    }
}

