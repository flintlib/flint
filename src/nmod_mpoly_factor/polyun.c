/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

/* assuming B is bivariate in x and y, put x on the outside, y inside */
void nmod_mpoly_get_polyu1n(
    n_polyun_t A,
    const nmod_mpoly_t B,
    slong varx,
    slong vary,
    const nmod_mpoly_ctx_t ctx)
{
    slong j, Ai;
    ulong Bexpx, Bexpy;
    slong Boffx, Bshiftx, Boffy, Bshifty;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    slong NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    
    mpoly_gen_offset_shift_sp(&Boffx, &Bshiftx, varx, B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&Boffy, &Bshifty, vary, B->bits, ctx->minfo);

    Ai = -1;
    for (j = 0; j < B->length; j++)
    {
        Bexpx = ((B->exps + NB*j)[Boffx] >> Bshiftx) & mask;
        Bexpy = ((B->exps + NB*j)[Boffy] >> Bshifty) & mask;

        if (Ai < 0 || A->exps[Ai] != Bexpx)
        {
            Ai++;
            n_polyun_fit_length(A, Ai + 1);
            A->exps[Ai] = Bexpx;
            n_poly_zero(A->coeffs + Ai);
        }

        n_poly_set_coeff(A->coeffs + Ai, Bexpy, B->coeffs[j]);
        if (n_poly_is_zero(A->coeffs + Ai))
            Ai--;
    }

    A->length = Ai + 1;

    FLINT_ASSERT(n_polyun_mod_is_canonical(A, ctx->mod));
}

void nmod_mpoly_set_polyu1n(
    nmod_mpoly_t B,
    const n_polyun_t A,
    slong varx,
    slong vary,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong Boffx, Bshiftx, Boffy, Bshifty;
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);

    mpoly_gen_offset_shift_sp(&Boffx, &Bshiftx, varx, B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&Boffy, &Bshifty, vary, B->bits, ctx->minfo);

    B->length = 0;
    for (i = 0; i < A->length; i++)
    {
        for (j = A->coeffs[i].length - 1; j >= 0; j--)
        {
            if (A->coeffs[i].coeffs[j] == 0)
                continue;

            nmod_mpoly_fit_length(B, B->length + 1, ctx);
            mpoly_monomial_zero(B->exps + N*B->length, N);
            (B->exps + N*B->length)[Boffx] += A->exps[i] << Bshiftx;
            (B->exps + N*B->length)[Boffy] += j << Bshifty;
            B->coeffs[B->length] = A->coeffs[i].coeffs[j];
            B->length++;
        }
    }

    FLINT_ASSERT(nmod_mpoly_is_canonical(B, ctx));
}

