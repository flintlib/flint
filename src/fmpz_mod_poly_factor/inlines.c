/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define FMPZ_MOD_POLY_FACTOR_INLINES_C

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_factor_get_fmpz_mod_poly(fmpz_mod_poly_t z,
                 fmpz_mod_poly_factor_t fac, slong i, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_set(z, fac->poly + i, ctx);
}

