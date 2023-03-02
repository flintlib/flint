/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"
#include "fmpz_mod_mpoly.h"
#include "fmpz_mod_mpoly_factor.h"

int fmpz_mod_mpoly_gcd_cofactors(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mod_mpoly_is_zero(B, ctx))
        {
            fmpz_mod_mpoly_zero(G, ctx);
            fmpz_mod_mpoly_zero(Abar, ctx);
            fmpz_mod_mpoly_zero(Bbar, ctx);
            return 1;
        }
        fmpz_mod_mpoly_set(G, B, ctx);
        fmpz_mod_mpoly_zero(Abar, ctx);
        fmpz_mod_mpoly_one(Bbar, ctx);
        if (!fmpz_is_one(G->coeffs + 0))
        {
            _fmpz_mod_vec_scalar_mul_fmpz_mod(Bbar->coeffs, Bbar->coeffs,
                                     Bbar->length, G->coeffs + 0, ctx->ffinfo);
            _fmpz_mod_vec_scalar_div_fmpz_mod(G->coeffs, G->coeffs,
                                        G->length, G->coeffs + 0, ctx->ffinfo);
        }
        return 1;
    }

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        fmpz_mod_mpoly_set(G, A, ctx);
        fmpz_mod_mpoly_zero(Bbar, ctx);
        fmpz_mod_mpoly_one(Abar, ctx);
        if (!fmpz_is_one(G->coeffs + 0))
        {
            _fmpz_mod_vec_scalar_mul_fmpz_mod(Abar->coeffs, Abar->coeffs,
                                     Abar->length, G->coeffs + 0, ctx->ffinfo);
            _fmpz_mod_vec_scalar_div_fmpz_mod(G->coeffs, G->coeffs,
                                        G->length, G->coeffs + 0, ctx->ffinfo);
        }
        return 1;
    }

    return _fmpz_mod_mpoly_gcd_algo(G, Abar, Bbar, A, B, ctx, MPOLY_GCD_USE_ALL);
}
