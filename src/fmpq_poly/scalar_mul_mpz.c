/*
    Copyright (C) 2010 Sebastian Pancratz

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

void 
fmpq_poly_scalar_mul_mpz(fmpq_poly_t rop, const fmpq_poly_t op, const mpz_t c)
{
    fmpz_t f;

    fmpz_init_set_readonly(f, c);
    fmpq_poly_scalar_mul_fmpz(rop, op, f);
    fmpz_clear_readonly(f);
}

