/*
    Copyright (C) 2018-2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"

int fmpz_mpoly_gcd(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
            fmpz_mpoly_zero(G, ctx);
        else if (fmpz_sgn(B->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, B, ctx);
        else
            fmpz_mpoly_set(G, B, ctx);

        return 1;
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, A, ctx);
        else
            fmpz_mpoly_set(G, A, ctx);

        return 1;
    }

    return _fmpz_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_ALL);
}

int fmpz_mpoly_gcd_brown(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx) || fmpz_mpoly_is_zero(B, ctx))
        return fmpz_mpoly_gcd(G, A, B, ctx);

    return _fmpz_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_BROWN);
}

int fmpz_mpoly_gcd_cofactors(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            fmpz_mpoly_zero(Abar, ctx);
            fmpz_mpoly_zero(Bbar, ctx);
            return 1;
        }
        fmpz_mpoly_set(G, B, ctx);
        fmpz_mpoly_zero(Abar, ctx);
        fmpz_mpoly_one(Bbar, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, G, ctx);
            fmpz_mpoly_neg(Bbar, Bbar, ctx);
        }
        return 1;
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        fmpz_mpoly_set(G, A, ctx);
        fmpz_mpoly_zero(Bbar, ctx);
        fmpz_mpoly_one(Abar, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, G, ctx);
            fmpz_mpoly_neg(Abar, Abar, ctx);
        }
        return 1;
    }

    return _fmpz_mpoly_gcd_algo(G, Abar, Bbar, A, B, ctx, MPOLY_GCD_USE_ALL);
}

int fmpz_mpoly_gcd_hensel(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx) || fmpz_mpoly_is_zero(B, ctx))
        return fmpz_mpoly_gcd(G, A, B, ctx);

    return _fmpz_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_HENSEL);
}

int fmpz_mpoly_gcd_subresultant(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx) || fmpz_mpoly_is_zero(B, ctx))
        return fmpz_mpoly_gcd(G, A, B, ctx);

    return _fmpz_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_PRS);
}

int fmpz_mpoly_gcd_zippel(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx) || fmpz_mpoly_is_zero(B, ctx))
        return fmpz_mpoly_gcd(G, A, B, ctx);

    return _fmpz_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_ZIPPEL);
}

int fmpz_mpoly_gcd_zippel2(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx) || fmpz_mpoly_is_zero(B, ctx))
        return fmpz_mpoly_gcd(G, A, B, ctx);

    return _fmpz_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_ZIPPEL2);
}
