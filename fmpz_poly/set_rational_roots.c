/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpq.h"
#include "fmpz_poly.h"

void fmpz_poly_set_rational_roots(fmpz_poly_t p, fmpq * vec, slong len)
{
    fmpz_poly_t q;
    slong i;

    fmpz_poly_init(q);
    fmpz_poly_one(p);
    for (i = 0; i < len; i++)
    {
        fmpz_poly_set_coeff_fmpz(q, 1, fmpq_denref(vec + i));
        fmpz_poly_set_coeff_fmpz(q, 0, fmpq_numref(vec + i));
        fmpz_neg(q->coeffs, q->coeffs);
        fmpz_poly_mul(p, p, q);
    }
    fmpz_poly_clear(q);
}
