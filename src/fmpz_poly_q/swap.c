/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_q.h"

void fmpz_poly_q_swap(fmpz_poly_q_t op1, fmpz_poly_q_t op2)
{
    if (op1 != op2)
    {
        fmpz_poly_struct *t;

        t        = op1->num;
        op1->num = op2->num;
        op2->num = t;

        t        = op1->den;
        op1->den = op2->den;
        op2->den = t;
    }
}
