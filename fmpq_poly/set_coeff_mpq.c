/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void fmpq_poly_set_coeff_mpq(fmpq_poly_t poly, slong n, const mpq_t x)
{
    fmpq_t f;

    fmpq_init_set_readonly(f, x);
    fmpq_poly_set_coeff_fmpq(poly, n, f);
    fmpq_clear_readonly(f);
}

