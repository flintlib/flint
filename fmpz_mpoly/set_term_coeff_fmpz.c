/*
    Copyright (C) 2008, 2009, 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_set_term_coeff_fmpz(fmpz_mpoly_t A,
                           slong i, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_set_term_coeff_fmpz");
    }

    fmpz_set(A->coeffs + i, c);
}

void fmpz_mpoly_set_term_coeff_ui(fmpz_mpoly_t A,
                                  slong i, ulong c, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_set_term_coeff_ui");
    }

    fmpz_set_ui(A->coeffs + i, c);
}

void fmpz_mpoly_set_term_coeff_si(fmpz_mpoly_t A,
                                  slong i, slong c, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_set_term_coeff_si");
    }

    fmpz_set_si(A->coeffs + i, c);
}
