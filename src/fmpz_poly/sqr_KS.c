/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_sqr_KS(fmpz *rop, const fmpz *op, slong len)
{
    _fmpz_poly_mul_KS(rop, op, len, op, len);
}

void fmpz_poly_sqr_KS(fmpz_poly_t rop, const fmpz_poly_t op)
{
    slong len;

    if (op->length == 0)
    {
        fmpz_poly_zero(rop);
        return;
    }

    len = 2 * op->length - 1;
    fmpz_poly_fit_length(rop, len);
    _fmpz_poly_sqr_KS(rop->coeffs, op->coeffs, op->length);
    _fmpz_poly_set_length(rop, len);
}
