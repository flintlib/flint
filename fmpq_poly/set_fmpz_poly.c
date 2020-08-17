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
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void fmpq_poly_set_fmpz_poly(fmpq_poly_t rop, const fmpz_poly_t op)
{
    if (fmpz_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        fmpq_poly_fit_length(rop, fmpz_poly_length(op));
        _fmpq_poly_set_length(rop, fmpz_poly_length(op));
        _fmpz_vec_set(rop->coeffs, op->coeffs, rop->length);
        fmpz_one(rop->den);
    }
}

