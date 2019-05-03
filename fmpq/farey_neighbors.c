/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

/*
    Find the fractions directly below and above a in the Farey sequence of
    order a_den:

     left_num     a_num = left_num + right_num     right_num
    ---------- < ------------------------------ < -----------
     left_den     a_den = left_den + right_den     right_den

    This function will return 0 if a is not canonical or has denominator 1.
    Otherwise it returns 1;
*/
int fmpq_farey_neighbors(fmpq_t left, fmpq_t right, const fmpq_t a)
{
    int success;

    FLINT_ASSERT(left != right);
    FLINT_ASSERT(a != right);
    FLINT_ASSERT(a != left);

    success = fmpz_invmod(fmpq_denref(left), fmpq_numref(a), fmpq_denref(a));
    if (!success || fmpz_is_zero(fmpq_denref(left)))
    {
        return 0;
    }

    fmpz_mul(fmpq_numref(left), fmpq_numref(a), fmpq_denref(left));
    fmpz_sub_ui(fmpq_numref(left), fmpq_numref(left), 1);
    fmpz_divexact(fmpq_numref(left), fmpq_numref(left), fmpq_denref(a));

    fmpz_sub(fmpq_numref(right), fmpq_numref(a), fmpq_numref(left));
    fmpz_sub(fmpq_denref(right), fmpq_denref(a), fmpq_denref(left));

    return 1;
}
