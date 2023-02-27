/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"

void fmpq_poly_set_coeff_mpz(fmpq_poly_t poly, slong n, const mpz_t x)
{
    fmpz_t f;

    fmpz_init_set_readonly(f, x);
    fmpq_poly_set_coeff_fmpz(poly, n, f);
    fmpz_clear_readonly(f);
}

