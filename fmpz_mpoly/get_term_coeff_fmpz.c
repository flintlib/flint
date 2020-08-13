/*
    Copyright (C) 2008, 2009, 2016 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_get_term_coeff_fmpz(fmpz_t c, const fmpz_mpoly_t A,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "index out of range in fmpz_mpoly_get_term_coeff_fmpz");
    }

    fmpz_set(c, A->coeffs + i);
}

ulong fmpz_mpoly_get_term_coeff_ui(const fmpz_mpoly_t A,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "index out of range in fmpz_mpoly_get_term_coeff_ui");
    }

    return fmpz_get_ui(A->coeffs + i);
}

slong fmpz_mpoly_get_term_coeff_si(const fmpz_mpoly_t A,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "index out of range in fmpz_mpoly_get_term_coeff_si");
    }

    return fmpz_get_si(A->coeffs + i);
}
