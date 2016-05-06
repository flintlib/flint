/*
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_factor(fmpz_mod_poly_factor_t res,
                     const fmpz_mod_poly_t f)
{
    slong n = fmpz_mod_poly_degree(f);
    mp_bitcnt_t bits = fmpz_bits(&f->p);

    if (5 * bits + n > 75)
        fmpz_mod_poly_factor_kaltofen_shoup(res, f);
    else
        fmpz_mod_poly_factor_cantor_zassenhaus(res, f);
}
