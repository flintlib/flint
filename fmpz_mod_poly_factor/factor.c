/*
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_factor(fmpz_mod_poly_factor_t res,
                             const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    slong n = fmpz_mod_poly_degree(f, ctx);
    flint_bitcnt_t bits = fmpz_bits(fmpz_mod_ctx_modulus(ctx));

    if (5 * bits + n > 75)
        fmpz_mod_poly_factor_kaltofen_shoup(res, f, ctx);
    else
        fmpz_mod_poly_factor_cantor_zassenhaus(res, f, ctx);
}
