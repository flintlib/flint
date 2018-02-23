/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


slong fmpz_mpoly_get_term_si_ui(const fmpz_mpoly_t poly,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t newc;
    slong ret;

    fmpz_init(newc);
    fmpz_mpoly_get_term_fmpz_ui(newc, poly, exp, ctx);

    ret = fmpz_get_si(newc);
    fmpz_clear(newc);
    return ret;
}
