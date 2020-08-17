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

void fmpq_poly_scalar_mul_mpq(fmpq_poly_t rop, const fmpq_poly_t op, const mpq_t c)
{
    fmpq_t f;

    fmpq_init_set_readonly(f, c);
    fmpq_poly_scalar_mul_fmpq(rop, op, f);
    fmpq_clear_readonly(f);
}

