/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_q.h"

void fmpz_poly_q_scalar_mul_mpz(fmpz_poly_q_t rop, 
                                const fmpz_poly_q_t op, const mpz_t x)
{
    fmpz_t y;

    fmpz_init(y);
    fmpz_set_mpz(y, x);

    fmpz_poly_scalar_mul_fmpz(rop->num, op->num, y);
    fmpz_poly_set(rop->den, op->den);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(y);
}
