/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly_q.h"

void
fmpz_poly_q_scalar_div_fmpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpq_t x)
{
    fmpz_t num, den;

    if (fmpz_sgn(fmpq_numref(x)) == 0)
        flint_throw(FLINT_ERROR, "Division by zero in %s\n", __func__);

    fmpz_init(num);
    fmpz_init(den);
    fmpz_set(num, fmpq_numref(x));
    fmpz_set(den, fmpq_denref(x));

    fmpz_poly_scalar_mul_fmpz(rop->num, op->num, den);
    fmpz_poly_scalar_mul_fmpz(rop->den, op->den, num);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(num);
    fmpz_clear(den);
}
