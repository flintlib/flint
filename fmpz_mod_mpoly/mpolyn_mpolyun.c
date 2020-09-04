/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

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
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
            A->coeffs = (fmpz_mod_poly_struct *) flint_malloc(new_alloc*sizeof(fmpz_mod_poly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
            A->coeffs = (fmpz_mod_poly_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz_mod_poly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_poly_init(A->coeffs + i, ctx->ffinfo);
        }
        A->alloc = new_alloc;
    }
}




void fmpz_mod_mpolyun_init(
    fmpz_mod_mpolyun_t A,
    flint_bitcnt_t bits,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fmpz_mod_mpolyun_clear(
    fmpz_mod_mpolyun_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mod_mpolyn_clear(A->coeffs + i, ctx);
    if (A->coeffs)
        flint_free(A->coeffs);
    if (A->exps)
        flint_free(A->exps);
}

/*
    get the leading coeff in x_0,...,x_var
    A is in R[x_0, ... x_(var-1)][x_var]
*/

fmpz * fmpz_mod_mpolyn_leadcoeff_last_ref(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_poly_struct * leadpoly;
    FLINT_ASSERT(A->length > 0);
    leadpoly = A->coeffs + 0;
    FLINT_ASSERT(leadpoly->length > 0);
    return leadpoly->coeffs + leadpoly->length - 1;
}

fmpz * fmpz_mod_mpolyun_leadcoeff_last_ref(
    const fmpz_mod_mpolyun_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return fmpz_mod_mpolyn_leadcoeff_last_ref(A->coeffs + 0, ctx);
}

fmpz_mod_poly_struct * fmpz_mod_mpolyn_leadcoeff_ref(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

fmpz_mod_poly_struct * fmpz_mod_mpolyun_leadcoeff_ref(
    fmpz_mod_mpolyun_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return fmpz_mod_mpolyn_leadcoeff_ref(A->coeffs + 0, ctx);
}

void fmpz_mod_mpolyun_swap(fmpz_mod_mpolyun_t A, fmpz_mod_mpolyun_t B)
{
   fmpz_mod_mpolyun_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mod_mpolyun_zero(
    fmpz_mod_mpolyun_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    A->length = 0;
}

void fmpz_mod_mpolyun_print_pretty(
    const fmpz_mod_mpolyun_t poly,
    const char ** x,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        FLINT_ASSERT((poly->coeffs + i)->bits == poly->bits);
        fmpz_mod_mpolyn_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fmpz_mod_mpolyun_fit_length(
    fmpz_mod_mpolyun_t A,
    slong length,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_mpolyn_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mod_mpolyn_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_mpolyn_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpz_mod_mpolyn_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_mpolyn_init(A->coeffs + i, A->bits, ctx);
        }
        A->alloc = new_alloc;
    }
}


void fmpz_mod_mpolyn_content_poly(
    fmpz_mod_poly_t a,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_t t;

    fmpz_mod_poly_init(t, ctx->ffinfo);

    fmpz_mod_poly_zero(a, ctx->ffinfo);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_gcd(t, a, B->coeffs + i, ctx->ffinfo);
        fmpz_mod_poly_swap(t, a, ctx->ffinfo);
        if (fmpz_mod_poly_degree(a, ctx->ffinfo) == 0)
            break;
    }

    fmpz_mod_poly_clear(t, ctx->ffinfo);
}

void fmpz_mod_mpolyun_content_last(
    fmpz_mod_poly_t a,
    const fmpz_mod_mpolyun_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_poly_t t;

    fmpz_mod_poly_init(t, ctx->ffinfo);

    fmpz_mod_poly_zero(a, ctx->ffinfo);
    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < (B->coeffs + i)->length; j++)
        {
            fmpz_mod_poly_gcd(t, a, (B->coeffs + i)->coeffs + j, ctx->ffinfo);
            fmpz_mod_poly_swap(t, a, ctx->ffinfo);
            if (fmpz_mod_poly_degree(a, ctx->ffinfo) == 0)
                break;
        }
    }

    fmpz_mod_poly_clear(t, ctx->ffinfo);
}


void fmpz_mod_mpolyn_divexact_poly(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_poly_t b,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_t r, t;

    fmpz_mod_poly_init(r, ctx->ffinfo);
    fmpz_mod_poly_init(t, ctx->ffinfo);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_divrem(t, r, A->coeffs + i, b, ctx->ffinfo);
        FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx->ffinfo));
        FLINT_ASSERT(!fmpz_mod_poly_is_zero(t, ctx->ffinfo));
        fmpz_mod_poly_swap(t, A->coeffs + i, ctx->ffinfo);
    }

    fmpz_mod_poly_clear(r, ctx->ffinfo);
    fmpz_mod_poly_clear(t, ctx->ffinfo);
}

void fmpz_mod_mpolyun_divexact_last(
    fmpz_mod_mpolyun_t A,
    const fmpz_mod_poly_t b,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_poly_t r, t;

    fmpz_mod_poly_init(r, ctx->ffinfo);
    fmpz_mod_poly_init(t, ctx->ffinfo);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_struct * Ac = (A->coeffs + i)->coeffs;
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fmpz_mod_poly_divrem(t, r, Ac + j, b, ctx->ffinfo);
            FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx->ffinfo));
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(t, ctx->ffinfo));
            fmpz_mod_poly_swap(t, Ac + j, ctx->ffinfo);
        }
    }
    fmpz_mod_poly_clear(r, ctx->ffinfo);
    fmpz_mod_poly_clear(t, ctx->ffinfo);
}


void fmpz_mod_mpolyn_mul_poly(
    fmpz_mod_mpolyn_t A,
    fmpz_mod_poly_t b,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_t t;

    fmpz_mod_poly_init(t, ctx->ffinfo);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_mul(t, A->coeffs + i, b, ctx->ffinfo);
        fmpz_mod_poly_swap(t, A->coeffs + i, ctx->ffinfo);
    }

    fmpz_mod_poly_clear(t, ctx->ffinfo);
}

void fmpz_mod_mpolyun_mul_last(
    fmpz_mod_mpolyun_t A,
    fmpz_mod_poly_t b,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_poly_t t;

    fmpz_mod_poly_init(t, ctx->ffinfo);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fmpz_mod_poly_mul(t, (A->coeffs + i)->coeffs + j, b, ctx->ffinfo);
            fmpz_mod_poly_swap(t, (A->coeffs + i)->coeffs + j, ctx->ffinfo);
        }
    }

    fmpz_mod_poly_clear(t, ctx->ffinfo);
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


slong fmpz_mod_mpolyun_lastdeg(
    const fmpz_mod_mpolyun_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong deg = -WORD(1);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            slong newdeg = fmpz_mod_poly_degree((A->coeffs + i)->coeffs + j, ctx->ffinfo);
            deg = FLINT_MAX(deg, newdeg);
        }
    }

    return deg;
}


void fmpz_mod_mpolyn_one(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N;

    fmpz_mod_mpolyn_fit_length(A, 1, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    fmpz_mod_poly_set_ui(Acoeff + 0, 1, ctx->ffinfo);
    mpoly_monomial_zero(Aexp + N*0, N);

    A->length = 1;
}

void fmpz_mod_mpolyun_one(
    fmpz_mod_mpolyun_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpolyun_fit_length(A, 1, ctx);
    fmpz_mod_mpolyn_one(A->coeffs + 0, ctx);
    A->exps[0] = 0;
    A->length = 1;
}





void fmpz_mod_mpolyn_scalar_mul_fmpz_mod(
    fmpz_mod_mpolyn_t A,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_scalar_mul_fmpz(A->coeffs + i, A->coeffs + i, c, ctx->ffinfo);
    }
}

void fmpz_mod_mpolyun_scalar_mul_fmpz_mod(
    fmpz_mod_mpolyun_t A,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    FLINT_ASSERT(!fmpz_is_zero(c));
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mpolyn_scalar_mul_fmpz_mod(A->coeffs + i, c, ctx);
    }
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
        {
            return 0;
        }
        if (!fmpz_mod_poly_equal(A->coeffs + i, B->coeffs + i, ctx->ffinfo))
        {
            return 0;
        }
    }
    return 1;
}

int fmpz_mod_mpolyun_equal(
    const fmpz_mod_mpolyun_t A,
    const fmpz_mod_mpolyun_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(A->bits == B->bits);

    if (A->length != B->length)
    {
        return 0;
    }
    for (i = 0; i < A->length; i++)
    {
        if (A->exps[i] != B->exps[i])
        {
            return 0;
        }
        if (!fmpz_mod_mpolyn_equal(A->coeffs + i, B->coeffs + i, ctx))
        {
            return 0;
        }
    }
    return 1;
}

