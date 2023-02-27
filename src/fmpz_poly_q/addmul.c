/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_q.h"

void fmpz_poly_q_addmul(fmpz_poly_q_t rop, 
                        const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
{
    fmpz_poly_q_t temp;

    fmpz_poly_q_init(temp);
    fmpz_poly_q_mul(temp, op1, op2);
    fmpz_poly_q_add(rop, rop, temp);
    fmpz_poly_q_clear(temp);
}
