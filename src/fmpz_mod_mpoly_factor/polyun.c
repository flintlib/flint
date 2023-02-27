/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"
#include "mpn_extras.h"
#include "fmpz_vec.h"


int fmpz_mod_polyun_is_canonical(
    const fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_mod_poly_is_canonical(A->coeffs + i, ctx) ||
            fmpz_mod_poly_is_zero(A->coeffs + i, ctx))
        {
            return 0;
        }
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}

void fmpz_mod_polyun_clear(
    fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    for (i = 0; i < A->alloc; i++)
        fmpz_mod_poly_clear(A->coeffs + i, ctx);

    flint_free(A->coeffs);
    flint_free(A->exps);
}

void fmpz_mod_polyun_realloc(
    fmpz_mod_polyun_t A,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    slong i, old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    A->exps = FLINT_ARRAY_REALLOC(A->exps, new_alloc, ulong);
    A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, new_alloc, fmpz_mod_poly_struct);

    for (i = old_alloc; i < new_alloc; i++)
        fmpz_mod_poly_init(A->coeffs + i, ctx);

    A->alloc = new_alloc;
}


void fmpz_mod_polyu2n_print_pretty(
    const fmpz_mod_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, varlast, ctx);
        flint_printf(")*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 1, 2),
            var1, extract_exp(A->exps[i], 0, 2));
    }

    if (first)
        flint_printf("0");
}

void fmpz_mod_polyu1n_print_pretty(
    const fmpz_mod_polyun_t A,
    const char * var0,
    const char * varlast,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, varlast, ctx);
        flint_printf(")*%s^%wu", var0, A->exps[i]);
    }

    if (first)
        flint_printf("0");
}


int fmpz_mod_polyun_equal(
    fmpz_mod_polyun_t A,
    const fmpz_mod_polyun_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (A->exps[i] != B->exps[i])
            return 0;
        if (!fmpz_mod_poly_equal(A->coeffs + i, B->coeffs + i, ctx))
            return 0;
    }
    return 1;
}

void fmpz_mod_polyun_set(
    fmpz_mod_polyun_t A,
    const fmpz_mod_polyun_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_polyun_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fmpz_mod_poly_set(A->coeffs + i, B->coeffs + i, ctx);
    }
    A->length = B->length;
}


void fmpz_mod_polyu3n_print_pretty(
    const fmpz_mod_polyun_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char * varlast,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, varlast, ctx);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}


void fmpz_mod_polyun_one(fmpz_mod_polyun_t A, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_polyun_fit_length(A, 1, ctx);
    fmpz_mod_poly_one(A->coeffs + 0, ctx);
    A->exps[0] = 0;
    A->length = 1;
}


void fmpz_mod_mpoly_get_polyu1n(
    fmpz_mod_polyun_t A,
    const fmpz_mod_mpoly_t B,
    slong varx,
    slong vary,
    const fmpz_mod_mpoly_ctx_t ctx)
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
            fmpz_mod_polyun_fit_length(A, Ai + 1, ctx->ffinfo);
            A->exps[Ai] = Bexpx;
            fmpz_mod_poly_zero(A->coeffs + Ai, ctx->ffinfo);
        }

        fmpz_mod_poly_set_coeff_fmpz(A->coeffs + Ai, Bexpy, B->coeffs + j, ctx->ffinfo);
        if (fmpz_mod_poly_is_zero(A->coeffs + Ai, ctx->ffinfo))
            Ai--;
    }

    A->length = Ai + 1;

    FLINT_ASSERT(fmpz_mod_polyun_is_canonical(A, ctx->ffinfo));
}



void fmpz_mod_mpoly_set_polyu1n(
    fmpz_mod_mpoly_t B,
    const fmpz_mod_polyun_t A,
    slong varx,
    slong vary,
    const fmpz_mod_mpoly_ctx_t ctx)
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
            if (fmpz_is_zero(A->coeffs[i].coeffs + j))
                continue;

            fmpz_mod_mpoly_fit_length(B, B->length + 1, ctx);
            mpoly_monomial_zero(B->exps + N*B->length, N);
            (B->exps + N*B->length)[Boffx] += A->exps[i] << Bshiftx;
            (B->exps + N*B->length)[Boffy] += j << Bshifty;
            fmpz_set(B->coeffs + B->length, A->coeffs[i].coeffs + j);
            B->length++;
        }
    }

    FLINT_ASSERT(fmpz_mod_mpoly_is_canonical(B, ctx));
}

