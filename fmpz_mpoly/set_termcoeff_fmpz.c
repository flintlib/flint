/*
    Copyright (C) 2008, 2009, 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_set_termcoeff_fmpz(fmpz_mpoly_t poly,
                           slong n, const fmpz_t x, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) n >= (ulong) poly->length)
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_termcoeff_fmpz");

    fmpz_set(poly->coeffs + n, x);
}

void fmpz_mpoly_set_termcoeff_ui(fmpz_mpoly_t poly,
                                  slong n, ulong x, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) n >= (ulong) poly->length)
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_termcoeff_ui");

    fmpz_set_ui(poly->coeffs + n, x);
}

void fmpz_mpoly_set_termcoeff_si(fmpz_mpoly_t poly,
                                  slong n, slong x, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) n >= (ulong) poly->length)
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_termcoeff_si");

    fmpz_set_si(poly->coeffs + n, x);
}
