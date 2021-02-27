/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_set_term_coeff_fmpz(fmpz_mod_mpoly_t A,
                      slong i, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (i >= (ulong) A->length)
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_set_term_coeff_fmpz: index is out of range");

    fmpz_mod_set_fmpz(A->coeffs + i, c, ctx->ffinfo);
}

void fmpz_mod_mpoly_set_term_coeff_ui(fmpz_mod_mpoly_t A,
                              slong i, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (i >= (ulong) A->length)
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_set_term_coeff_ui: index is out of range");

    fmpz_mod_set_ui(A->coeffs + i, c, ctx->ffinfo);
}

void fmpz_mod_mpoly_set_term_coeff_si(fmpz_mod_mpoly_t A,
                              slong i, slong c, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (i >= (ulong) A->length)
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_set_term_coeff_si: index is out of range");

    fmpz_mod_set_si(A->coeffs + i, c, ctx->ffinfo);
}
