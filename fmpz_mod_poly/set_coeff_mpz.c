/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_set_coeff_mpz(fmpz_mod_poly_t poly, slong n, const mpz_t x,
                                                    const fmpz_mod_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_readonly(t, x);
    fmpz_mod_poly_set_coeff_fmpz(poly, n, t, ctx);
    fmpz_clear_readonly(t);
}
