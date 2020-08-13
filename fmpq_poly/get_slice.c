/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void fmpq_poly_get_slice(fmpq_poly_t rop, const fmpq_poly_t op, slong i, slong j)
{
    i = FLINT_MAX(i, 0);
    j = FLINT_MIN(j, op->length);

    if (i < j)
    {
        slong k;

        if (rop == op)
        {
            for (k = 0; k < i; k++)
                fmpz_zero(rop->coeffs + k);
            for (k = j; k < rop->length; k++)
                fmpz_zero(rop->coeffs + k);
            fmpq_poly_canonicalise(rop);
        }
        else
        {
            fmpq_poly_fit_length(rop, j);
            _fmpq_poly_set_length(rop, j);

            _fmpz_vec_set(rop->coeffs + i, op->coeffs + i, j - i);
            fmpz_set(rop->den, op->den);
            fmpq_poly_canonicalise(rop);
        }
    }
    else
    {
        fmpq_poly_zero(rop);
    }
}

