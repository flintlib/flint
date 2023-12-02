/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly_q.h"

void fmpz_poly_q_inv(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    if (fmpz_poly_is_zero(op->num))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_q_inv). Zero is not invertible.\n");
    }

    if (rop == op)
    {
        fmpz_poly_swap(rop->num, rop->den);
        if (fmpz_sgn(fmpz_poly_lead(rop->den)) < 0)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_poly_neg(rop->den, rop->den);
        }
    }
    else
    {
        if (fmpz_sgn(fmpz_poly_lead(op->num)) > 0)
        {
            fmpz_poly_set(rop->num, op->den);
            fmpz_poly_set(rop->den, op->num);
        }
        else
        {
            fmpz_poly_neg(rop->num, op->den);
            fmpz_poly_neg(rop->den, op->num);
        }
    }
}

