/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_q.h"

void fmpz_poly_q_pow(fmpz_poly_q_t rop, 
                     const fmpz_poly_q_t op, ulong exp)
{
    if (exp == 0)
    {
        fmpz_poly_q_one(rop);
    }
    else
    {
        fmpz_poly_pow(rop->num, op->num, exp);
        fmpz_poly_pow(rop->den, op->den, exp);
    }
}
