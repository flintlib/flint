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
fmpz_poly_q_scalar_div_fmpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpz_t x)
{
    fmpz_t y;

    if (fmpz_sgn(x) == 0)
        flint_throw(FLINT_ERROR, "Division by zero in %s\n", __func__);

    fmpz_init(y);
    fmpz_set(y, x);

    fmpz_poly_set(rop->num, op->num);
    fmpz_poly_scalar_mul_fmpz(rop->den, op->den, y);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(y);
}
