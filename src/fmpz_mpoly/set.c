/*
    Copyright (C) 2008, 2009, 2016 William Hart
    Copyright (C) 2017, 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void _fmpz_mpoly_set(fmpz * poly1, ulong * exps1,
                     const fmpz * poly2, const ulong * exps2, slong n, slong N)
{
   slong i;

   if (poly1 != poly2)
   {
      for (i = 0; i < n; i++)
         fmpz_set(poly1 + i, poly2 + i);
   }

   if (exps1 != exps2)
   {
      for (i = 0; i < n*N; i++)
         exps1[i] = exps2[i];
   }
}

void fmpz_mpoly_set(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong N;
   N = mpoly_words_per_exp(poly2->bits, ctx->minfo);

   fmpz_mpoly_fit_length(poly1, poly2->length, ctx);
   fmpz_mpoly_fit_bits(poly1, poly2->bits, ctx);

   _fmpz_mpoly_set(poly1->coeffs, poly1->exps,
                   poly2->coeffs, poly2->exps, poly2->length, N);

   _fmpz_mpoly_set_length(poly1, poly2->length, ctx);
   poly1->bits = poly2->bits;
}

void _fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t poly,
                  const fmpz_t c, const fmpz * exp, const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t exp_bits;
    slong i, N, index;
    ulong * cmpmask;
    ulong * packed_exp;
    int exists;
    TMP_INIT;

    TMP_START;

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mpoly_fit_bits(poly, exp_bits, ctx);

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, poly->bits, ctx->minfo);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_set_monomial_ffmpz(packed_exp, exp, poly->bits, ctx->minfo);
    exists = mpoly_monomial_exists(&index, poly->exps,
                                  packed_exp, poly->length, N, cmpmask);

    if (!exists)
    {
        if (!fmpz_is_zero(c)) /* make new term only if coeff is nonzero*/
        {       

            fmpz_mpoly_fit_length(poly, poly->length + 1, ctx);

            for (i = poly->length; i >= index + 1; i--)
            {
                fmpz_set(poly->coeffs + i, poly->coeffs + i - 1);
                mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i - 1), N);
            }

            fmpz_set(poly->coeffs + index, c);
            mpoly_monomial_set(poly->exps + N*index, packed_exp, N);

            poly->length++; /* safe because length is increasing */
        }
    } else if (fmpz_is_zero(c)) /* zero coeff, remove term */
    {
        for (i = index; i < poly->length - 1; i++)
        {
            fmpz_set(poly->coeffs + i, poly->coeffs + i + 1);
            mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i + 1), N);
        }

        _fmpz_mpoly_set_length(poly, poly->length - 1, ctx);

    } else /* term with that monomial exists, coeff is nonzero */
    {
        fmpz_set(poly->coeffs + index, c);  
    }

   TMP_END; 
}


void fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t poly,
                const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * newexp;
    TMP_INIT;

    TMP_START;
    newexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
    {
        fmpz_init(newexp + i);
        fmpz_set(newexp + i, exp[i]);
    }

    _fmpz_mpoly_set_coeff_fmpz_fmpz(poly, c, newexp, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;
}

void fmpz_mpoly_set_coeff_fmpz_monomial(fmpz_mpoly_t A, const fmpz_t c,
                         const fmpz_mpoly_t M, const fmpz_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * texps;
    TMP_INIT;

    if (M->length != WORD(1))
    {
        flint_throw(FLINT_ERROR, "M not monomial in fmpz_mpoly_set_coeff_fmpz_monomial");
    }

    TMP_START;
    texps = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
        fmpz_init(texps + i);

    mpoly_get_monomial_ffmpz(texps, M->exps + 0, M->bits, ctx->minfo);
    _fmpz_mpoly_set_coeff_fmpz_fmpz(A, c, texps, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(texps + i);

    TMP_END;
    return;
}

void fmpz_mpoly_set_coeff_fmpz_ui(fmpz_mpoly_t poly,
                 const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * newexp;
    TMP_INIT;

    TMP_START;
    newexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
        fmpz_init_set_ui(newexp + i, exp[i]);

    _fmpz_mpoly_set_coeff_fmpz_fmpz(poly, c, newexp, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;
}

void fmpz_mpoly_set_coeff_si_fmpz(fmpz_mpoly_t poly,
                 const slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t newc;

    fmpz_init(newc);
    fmpz_set_si(newc, c);
    fmpz_mpoly_set_coeff_fmpz_fmpz(poly, newc, exp, ctx);
    fmpz_clear(newc);
}

void fmpz_mpoly_set_coeff_si_ui(fmpz_mpoly_t poly,
                 const slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t newc;

    fmpz_init(newc);
    fmpz_set_si(newc, c);
    fmpz_mpoly_set_coeff_fmpz_ui(poly, newc, exp, ctx);
    fmpz_clear(newc);
}

void fmpz_mpoly_set_coeff_ui_fmpz(fmpz_mpoly_t poly,
                const ulong c, fmpz * const *  exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t newc;

    fmpz_init(newc);
    fmpz_set_ui(newc, c);
    fmpz_mpoly_set_coeff_fmpz_fmpz(poly, newc, exp, ctx);
    fmpz_clear(newc);
}

void fmpz_mpoly_set_coeff_ui_ui(fmpz_mpoly_t poly,
                 const ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t newc;

    fmpz_init(newc);
    fmpz_set_ui(newc, c);
    fmpz_mpoly_set_coeff_fmpz_ui(poly, newc, exp, ctx);
    fmpz_clear(newc);
}

void fmpz_mpoly_set_si(fmpz_mpoly_t A, slong c, const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (c == 0)
    {
        _fmpz_mpoly_set_length(A, 0, ctx);
        return;
    }

    fmpz_mpoly_fit_length(A, 1, ctx);
    fmpz_set_si(A->coeffs + 0, c);
    mpoly_monomial_zero(A->exps + N*0, N);
    _fmpz_mpoly_set_length(A, 1, ctx);
}

void fmpz_mpoly_set_fmpz(fmpz_mpoly_t A,
                                    const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (fmpz_is_zero(c))
    {
        _fmpz_mpoly_set_length(A, 0, ctx);
        return;
    }

    fmpz_mpoly_fit_length(A, 1, ctx);
    fmpz_set(A->coeffs + 0, c);
    mpoly_monomial_zero(A->exps + N*0, N);
    _fmpz_mpoly_set_length(A, 1, ctx);
}

void fmpz_mpoly_set_ui(fmpz_mpoly_t A, ulong c, const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (c == 0)
    {
        _fmpz_mpoly_set_length(A, 0, ctx);
        return;
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    fmpz_mpoly_fit_length(A, 1, ctx);
    fmpz_set_ui(A->coeffs + 0, c);
    mpoly_monomial_zero(A->exps + N*0, N);
    _fmpz_mpoly_set_length(A, 1, ctx);
}

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

void fmpz_mpoly_set_term_exp_fmpz(fmpz_mpoly_t A, 
                       slong i, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t exp_bits;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_set_term_exp_fmpz");
    }

    exp_bits = mpoly_exp_bits_required_pfmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mpoly_fit_bits(A, exp_bits, ctx);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_set_monomial_pfmpz(A->exps + N*i, exp, A->bits, ctx->minfo);
}

void fmpz_mpoly_set_term_exp_ui(fmpz_mpoly_t A, 
                       slong i, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t exp_bits;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_set_term_exp_ui");
    }

    exp_bits = mpoly_exp_bits_required_ui(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mpoly_fit_bits(A, exp_bits, ctx);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_set_monomial_ui(A->exps + N*i, exp, A->bits, ctx->minfo);
}
