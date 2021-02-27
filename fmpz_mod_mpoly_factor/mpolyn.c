/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"

void fmpz_mod_mpolyn_init(
    fmpz_mod_mpolyn_t A,
    flint_bitcnt_t bits,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fmpz_mod_mpolyn_clear(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mod_poly_clear(A->coeffs + i, ctx->ffinfo);
    if (A->coeffs)
        flint_free(A->coeffs);
    if (A->exps)
        flint_free(A->exps);
}

int fmpz_mod_mpolyn_is_canonical(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    if (!mpoly_monomials_valid_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    if (mpoly_monomials_overflow_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    if (!mpoly_monomials_inorder_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_mod_poly_is_canonical(A->coeffs + i, ctx->ffinfo) ||
            fmpz_mod_poly_is_zero(A->coeffs + i, ctx->ffinfo))
        {
            return 0;
        }
    }

    return 1;
}

void fmpz_mod_mpolyn_print_pretty(
    const fmpz_mod_mpolyn_t A,
    const char ** x_in,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_poly_struct * coeff = A->coeffs;
    slong len = A->length;
    ulong * exp = A->exps;
    slong bits = A->bits;
    slong i, j, N;
    fmpz * exponents;
    char ** x = (char **) x_in;
    TMP_INIT;

    if (len == 0)
    {
        flint_printf("0");
        return;
    }

    N = mpoly_words_per_exp(bits, ctx->minfo);

    TMP_START;

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(ctx->minfo->nvars*sizeof(char *));
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(((FLINT_BITS+4)/3)*sizeof(char));
            flint_sprintf(x[i], "x%wd", i+1);
        }
    }

    exponents = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(exponents + i);
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            printf(" + ");
        }

        printf("(");
        fmpz_mod_poly_print_pretty(coeff + i, "v", ctx->ffinfo);
        printf(")");

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, ctx->minfo);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, WORD(1));

            if (cmp > 0)
            {
                printf("*%s^", x[j]);
                fmpz_print(exponents + j);
            } else if (cmp == 0)
            {
                printf("*%s", x[j]);
            }
        }
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;
}

void fmpz_mod_mpolyn_fit_length(
    fmpz_mod_mpolyn_t A,
    slong length,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length <= old_alloc)
        return;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    A->exps = FLINT_ARRAY_REALLOC(A->exps, new_alloc*N, ulong);
    A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, new_alloc, fmpz_mod_poly_struct);

    for (i = old_alloc; i < new_alloc; i++)
        fmpz_mod_poly_init(A->coeffs + i, ctx->ffinfo);

    A->alloc = new_alloc;
}


slong fmpz_mod_mpolyn_lastdeg(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    slong deg = -WORD(1);

    for (i = 0; i < A->length; i++)
    {
        slong newdeg = fmpz_mod_poly_degree(A->coeffs + i, ctx->ffinfo);
        deg = FLINT_MAX(deg, newdeg);
    }

    return deg;
}


void fmpz_mod_mpolyn_one(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    fmpz_mod_mpolyn_fit_length(A, 1, ctx);
    fmpz_mod_poly_set_ui(A->coeffs + 0, 1, ctx->ffinfo);
    mpoly_monomial_zero(A->exps + N*0, N);
    A->length = 1;
}


void fmpz_mod_mpolyn_scalar_mul_fmpz_mod(
    fmpz_mod_mpolyn_t A,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_scalar_mul_fmpz(A->coeffs + i, A->coeffs + i, c, ctx->ffinfo);
}


int fmpz_mod_mpolyn_equal(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong i;

    FLINT_ASSERT(A->bits == B->bits);

    if (A->length != B->length)
    {
        return 0;
    }
    for (i = 0; i < A->length; i++)
    {
        if (!mpoly_monomial_equal(A->exps + N*i, B->exps + N*i, N))
            return 0;

        if (!fmpz_mod_poly_equal(A->coeffs + i, B->coeffs + i, ctx->ffinfo))
            return 0;
    }
    return 1;
}


/* put the last variable of B back into A */
void fmpz_mod_mpoly_cvtfrom_mpolyn(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j, k;
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    ulong * genexp;
    TMP_INIT;

    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    TMP_START;

    genexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(genexp, var, B->bits, ctx->minfo);

    fmpz_mod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);

    k = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = B->coeffs[i].length - 1; j >= 0; j--)
        {
            fmpz * c = B->coeffs[i].coeffs + j;
            if (fmpz_is_zero(c))
                continue;

            _fmpz_mod_mpoly_fit_length(&A->coeffs, &A->coeffs_alloc,
                                       &A->exps, &A->exps_alloc, N, k + 1);
            fmpz_set(A->coeffs + k, c);
            mpoly_monomial_madd(A->exps + N*k, B->exps + N*i, j, genexp, N);                
            k++;
        }
    }

    A->length = k;
    TMP_END;
}


/* take the last variable of B out */
void fmpz_mod_mpoly_cvtto_mpolyn(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    ulong * oneexp;
    slong offset;
    slong shift;
    ulong mask;
    slong N;
    TMP_INIT;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    TMP_START;

    N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    mpoly_gen_monomial_offset_shift_sp(oneexp, &offset, &shift, var,
                                                          B->bits, ctx->minfo);

    k = 0;
    fmpz_mod_mpolyn_fit_length(A, k + 1, ctx);
    for (i = 0; i < B->length; i++)
    {
        ulong c = (B->exps[N*i + offset] >> shift) & mask;
        mpoly_monomial_msub(A->exps + N*k, B->exps + N*i, c, oneexp, N);

        if (k > 0 && mpoly_monomial_equal(A->exps + N*k, A->exps + N*(k - 1), N))
        {
            fmpz_mod_poly_set_coeff_fmpz(A->coeffs + k - 1, c, B->coeffs + i, ctx->ffinfo);
        }
        else
        {
            fmpz_mod_poly_zero(A->coeffs + k, ctx->ffinfo);
            fmpz_mod_poly_set_coeff_fmpz(A->coeffs + k, c, B->coeffs + i, ctx->ffinfo);
            k++;
            fmpz_mod_mpolyn_fit_length(A, k + 1, ctx);
        }
    }

    A->length = k;
    TMP_END;
}


void fmpz_mod_mpoly_to_mpolyn_perm_deflate(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t nctx,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong j, k, l;
    slong NA = mpoly_words_per_exp_sp(A->bits, nctx->minfo);
    slong NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong n = ctx->minfo->nvars;
    slong m = nctx->minfo->nvars;
    ulong * Bexps;
    slong * offs, * shifts;
    fmpz_mod_mpoly_t T;
    TMP_INIT;

    FLINT_ASSERT(m <= n);

    TMP_START;
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    offs   = (slong *) TMP_ALLOC(m*sizeof(slong));
    shifts = (slong *) TMP_ALLOC(m*sizeof(slong));
    for (k = 0; k < m; k++)
    {
        mpoly_gen_offset_shift_sp(offs + k, shifts + k, k, A->bits, nctx->minfo);
    }

    fmpz_mod_mpoly_init3(T, B->length, A->bits, nctx);
    T->length = B->length;
    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        fmpz_set(T->coeffs + j, B->coeffs + j);
        mpoly_monomial_zero(T->exps + NA*j, NA);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            (T->exps + NA*j)[offs[k]] += ((Bexps[l] - shift[l]) / stride[l]) << shifts[k];
        }
    }

    fmpz_mod_mpoly_sort_terms(T, nctx);

    fmpz_mod_mpoly_cvtto_mpolyn(A, T, nctx->minfo->nvars - 1, nctx);

    fmpz_mod_mpoly_clear(T, nctx);

    TMP_END;
}

void fmpz_mod_mpoly_from_mpolyn_perm_inflate(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mod_mpoly_ctx_t ctx,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t nctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong n = ctx->minfo->nvars;
    slong m = nctx->minfo->nvars;
    slong i, h, k, l;
    slong NA, NB;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    ulong * Bexps;
    ulong * Aexps, * tAexp, * tAgexp;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m <= n);
    TMP_START;

    Bexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, nctx->minfo);

    tAexp = (ulong *) TMP_ALLOC(NA*sizeof(ulong));
    tAgexp = (ulong *) TMP_ALLOC(NA*sizeof(ulong));
    mpoly_gen_monomial_sp(tAgexp, perm[m - 1], Abits, ctx->minfo);
    for (i = 0; i < NA; i++)
        tAgexp[i] *= stride[perm[m - 1]];

    fmpz_mod_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, nctx->minfo);
        FLINT_ASSERT(Bexps[m - 1] == 0);
        for (l = 0; l < n; l++)
        {
            Aexps[l] = shift[l];
        }
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            Aexps[l] += stride[l]*Bexps[k];
        }

        mpoly_set_monomial_ui(tAexp, Aexps, Abits, ctx->minfo);

        h = (B->coeffs + i)->length;
        _fmpz_mod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                                   &Aexp, &A->exps_alloc, NA, Alen + h);
        for (h--; h >= 0; h--)
        {
            fmpz * c = (B->coeffs + i)->coeffs + h;
            if (fmpz_is_zero(c))
                continue;
            mpoly_monomial_madd(Aexp + NA*Alen, tAexp, h, tAgexp, NA);
            fmpz_set(Acoeff + Alen, c);
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    _fmpz_mod_mpoly_set_length(A, Alen, ctx);

    fmpz_mod_mpoly_sort_terms(A, ctx);
    TMP_END;
}


void fmpz_mod_mpolyn_set(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    slong Blen = B->length;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    FLINT_ASSERT(A->bits == B->bits);

    fmpz_mod_mpolyn_fit_length(A, Blen, ctx);

    mpoly_copy_monomials(A->exps, B->exps, Blen, N);

    for (i = 0; i < Blen; i++)
        fmpz_mod_poly_set(A->coeffs + i, B->coeffs + i, ctx->ffinfo);

    A->length = Blen;
}


int fmpz_mod_mpolyn_is_nonzero_fmpz(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length != 1 || A->coeffs[0].length != 1)
        return 0;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    return mpoly_monomial_is_zero(A->exps + N*0, N);
}

