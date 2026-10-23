/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

/* set res to f divided by its leading coefficient, or to zero if f is zero */
void
fmpz_mod_mpoly_monic_part(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t f, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_zero(f, ctx))
        fmpz_mod_mpoly_zero(res, ctx);
    else
        fmpz_mod_mpoly_make_monic(res, f, ctx);
}
