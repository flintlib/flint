/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_get_term_coeff_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A,
                                       slong i, const fmpz_mod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    if (i >= A->length)
        flint_throw(FLINT_ERROR,
                     "fmpz_mod_mpoly_get_term_coeff_fmpz: index out of range");

    fmpz_set(c, A->coeffs + i);
}
