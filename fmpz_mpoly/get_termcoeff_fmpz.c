/*
    Copyright (C) 2008, 2009, 2016 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_get_termcoeff_fmpz(fmpz_t x, const fmpz_mpoly_t poly,
                                           slong n, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) n >= (ulong) poly->length)
        flint_throw(FLINT_ERROR, "index out of range in fmpz_mpoly_get_termcoeff_fmpz");

    fmpz_set(x, poly->coeffs + n);
}

ulong fmpz_mpoly_get_termcoeff_ui(const fmpz_mpoly_t poly,
                                           slong n, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) n >= (ulong) poly->length)
        flint_throw(FLINT_ERROR, "index out of range in fmpz_mpoly_get_termcoeff_ui");

    return fmpz_get_ui(poly->coeffs + n);
}

slong fmpz_mpoly_get_termcoeff_si(const fmpz_mpoly_t poly,
                                           slong n, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) n >= (ulong) poly->length)
        flint_throw(FLINT_ERROR, "index out of range in fmpz_mpoly_get_termcoeff_si");

    return fmpz_get_si(poly->coeffs + n);
}
