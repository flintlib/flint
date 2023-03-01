/*
    Copyright (C) 2008, 2009, 2016 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2017-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t A,
                                fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong index;
    index = mpoly_monomial_index_pfmpz(A->exps, A->bits, A->length,
                                                              exp, ctx->minfo);
    if (index < 0)
    {
        fmpz_zero(c);
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        fmpz_set(c, A->coeffs + index);
    }
}

void fmpz_mpoly_get_coeff_fmpz_monomial(fmpz_t c, const fmpz_mpoly_t A,
                         const fmpz_mpoly_t M, const fmpz_mpoly_ctx_t ctx)
{
    slong index;

    if (M->length != WORD(1))
    {
        flint_throw(FLINT_ERROR, "M not monomial in fmpz_mpoly_get_coeff_fmpz_monomial");
    }

    index = mpoly_monomial_index_monomial(A->exps, A->bits, A->length,
                                                 M->exps, M->bits, ctx->minfo);
    if (index < 0)
    {
        fmpz_zero(c);
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        fmpz_set(c, A->coeffs + index);
    }
}

void fmpz_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mpoly_t A,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong index;
    index = mpoly_monomial_index_ui(A->exps, A->bits, A->length,
                                                              exp, ctx->minfo);
    if (index < 0)
    {
        fmpz_zero(c);
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        fmpz_set(c, A->coeffs + index);
    }
}

slong fmpz_mpoly_get_coeff_si_fmpz(const fmpz_mpoly_t A,
                                fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong index;
    index = mpoly_monomial_index_pfmpz(A->exps, A->bits, A->length,
                                                              exp, ctx->minfo);
    if (index < 0)
    {
        return 0;
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        return fmpz_get_si(A->coeffs + index);
    }
}

slong fmpz_mpoly_get_coeff_si_ui(const fmpz_mpoly_t A,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong index;
    index = mpoly_monomial_index_ui(A->exps, A->bits, A->length,
                                                              exp, ctx->minfo);
    if (index < 0)
    {
        return 0;
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        return fmpz_get_si(A->coeffs + index);
    }
}

ulong fmpz_mpoly_get_coeff_ui_fmpz(const fmpz_mpoly_t A,
                                fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong index;
    index = mpoly_monomial_index_pfmpz(A->exps, A->bits, A->length,
                                                              exp, ctx->minfo);
    if (index < 0)
    {
        return 0;
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        return fmpz_get_ui(A->coeffs + index);
    }
}

ulong fmpz_mpoly_get_coeff_ui_ui(const fmpz_mpoly_t A,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong index;
    index = mpoly_monomial_index_ui(A->exps, A->bits, A->length,
                                                              exp, ctx->minfo);
    if (index < 0)
    {
        return 0;
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        return fmpz_get_ui(A->coeffs + index);
    }
}

void fmpz_mpoly_get_coeff_vars_ui(fmpz_mpoly_t C, const fmpz_mpoly_t A,
                         const slong * vars, const ulong * exps, slong length,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, N;
    slong offset, shift;
    slong maxoffset, minoffset;
    ulong * uexp;
    ulong * tmask, * texp;
    slong nvars = ctx->minfo->nvars;
    fmpz * Ccoeff;
    ulong * Cexp;
    slong Calloc;
    slong Clen;
    TMP_INIT;

    if (C == A)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init(T, ctx);
        fmpz_mpoly_get_coeff_vars_ui(T, A, vars, exps,length, ctx);
        fmpz_mpoly_swap(T, C, ctx);
        fmpz_mpoly_clear(T, ctx);
        return;
    }

    TMP_START;

    uexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    for (i = 0; i < nvars; i++)
    {
        uexp[i] = 0;
    }
    for (i = 0; i < length; i++)
    {
        uexp[vars[i]] = exps[i];
    }

    if (A->bits < mpoly_exp_bits_required_ui(uexp, ctx->minfo))
    {
        fmpz_mpoly_zero(C, ctx);
        goto cleanup;
    }

    fmpz_mpoly_fit_bits(C, A->bits, ctx);
    C->bits = A->bits;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    tmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(tmask, N);
    mpoly_set_monomial_ui(texp, uexp, A->bits, ctx->minfo);

    if (A->bits <= FLINT_BITS)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
        maxoffset = 0;
        minoffset = N;
        for (i = 0; i < length; i++)
        {
            mpoly_gen_offset_shift_sp(&offset, &shift, vars[i], A->bits, ctx->minfo);
            tmask[offset] |= mask << shift;
            maxoffset = FLINT_MAX(maxoffset, offset);
            minoffset = FLINT_MIN(minoffset, offset);
        }
        FLINT_ASSERT(minoffset < N);

        Ccoeff = C->coeffs;
        Cexp = C->exps;
        Calloc = C->alloc;
        Clen = 0;
        for (i = 0; i < A->length; i++)
        {
            for (j = minoffset; j <= maxoffset; j++)
            {
                if ((((A->exps + N*i)[j] ^ texp[j]) & tmask[j]) != UWORD(0))
                    goto continue_outer_sp;
            }
            _fmpz_mpoly_fit_length(&Ccoeff, &Cexp, &Calloc, Clen + 1, N);
            mpoly_monomial_sub(Cexp + N*Clen, A->exps + N*i, texp, N);
            fmpz_set(Ccoeff + Clen, A->coeffs + i);
            Clen++;
continue_outer_sp:;
        }

        C->coeffs = Ccoeff;
        C->exps = Cexp;
        C->alloc = Calloc;
        _fmpz_mpoly_set_length(C, Clen, ctx);
    }
    else
    {
        ulong wpf = A->bits/FLINT_BITS;
        maxoffset = 0;
        minoffset = N;
        for (i = 0; i < length; i++)
        {
            offset = mpoly_gen_offset_mp(vars[i], A->bits, ctx->minfo);
            minoffset = FLINT_MIN(minoffset, offset);
            maxoffset = FLINT_MAX(maxoffset, offset + wpf - 1);
            for (j = 0; j < wpf; j++)
                tmask[offset + j] = -UWORD(1);
        }
        FLINT_ASSERT(minoffset < N);

        Ccoeff = C->coeffs;
        Cexp = C->exps;
        Calloc = C->alloc;
        Clen = 0;
        for (i = 0; i < A->length; i++)
        {
            for (j = minoffset; j <= maxoffset; j++)
            {
                if ((((A->exps + N*i)[j] ^ texp[j]) & tmask[j]) != UWORD(0))
                    goto continue_outer_mp;
            }
            _fmpz_mpoly_fit_length(&Ccoeff, &Cexp, &Calloc, Clen + 1, N);
            mpoly_monomial_sub_mp(Cexp + N*Clen, A->exps + N*i, texp, N);
            fmpz_set(Ccoeff + Clen, A->coeffs + i);
            Clen++;
continue_outer_mp:;
        }

        C->coeffs = Ccoeff;
        C->exps = Cexp;
        C->alloc = Calloc;
        _fmpz_mpoly_set_length(C, Clen, ctx);
    }

cleanup:

    TMP_END;
    return;
}

void fmpz_mpoly_get_fmpz(fmpz_t c, const fmpz_mpoly_t A, 
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length > WORD(1))
        flint_throw(FLINT_ERROR, "Nonconstant polynomial in fmpz_mpoly_get_fmpz");

    if (A->length == WORD(0))
    {
        fmpz_zero(c);
        return;
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    if (!mpoly_monomial_is_zero(A->exps + N*0, N))
        flint_throw(FLINT_ERROR, "Nonconstant monomial in fmpz_mpoly_get_fmpz");

    fmpz_set(c, A->coeffs + 0);
}

mpz_srcptr _fmpz_mpoly_get_mpz_signed_uiuiui(ulong * sm, fmpz x, mpz_ptr t)
{
    mpz_ptr p;
    slong i, abs_size;
    ulong s;

    if (!COEFF_IS_MPZ(x))
    {
        sm[0] = x;
        sm[1] = FLINT_SIGN_EXT(x);
        sm[2] = FLINT_SIGN_EXT(x);
    }
    else
    {
        p = COEFF_TO_PTR(x);

        sm[0] = 0;
        sm[1] = 0;
        sm[2] = 0;

        s = FLINT_SIGN_EXT(p->_mp_size);
        abs_size = FLINT_ABS(p->_mp_size);

        if (abs_size > 3 || (abs_size == 3 && p->_mp_d[2] >= COEFF_MAX))
            return p;

        for (i = 0; i < abs_size; i++)
            sm[i] = p->_mp_d[i];

        sub_dddmmmsss(sm[2], sm[1], sm[0], s^sm[2], s^sm[1], s^sm[0], s, s, s);
    }

    mpz_set_ui(t, 0);
    return t;
}

int fmpz_mpoly_is_fmpz_poly(
    const fmpz_mpoly_t A,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    return mpoly_is_poly(A->exps, A->length, A->bits, var, ctx->minfo);
}

int fmpz_mpoly_get_fmpz_poly(
    fmpz_poly_t A,
    const fmpz_mpoly_t B,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    slong Blen = B->length;
    const fmpz * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong i, N = mpoly_words_per_exp(Bbits, ctx->minfo);
    ulong k;

    fmpz_poly_zero(A);

    if (B->length < 1)
        return 1;

    if (Bbits <= FLINT_BITS)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - Bbits);
        slong off, shift;

        mpoly_gen_offset_shift_sp(&off, &shift, var, Bbits, ctx->minfo);

        for (i = 0; i < Blen; i++)
        {
            k = (Bexps[N*i + off] >> shift) & mask;
            fmpz_poly_set_coeff_fmpz(A, k, Bcoeffs + i);
        }
        return 1;
    }
    else
    {
        slong j, off;
        ulong check, wpf = Bbits/FLINT_BITS;

        off = mpoly_gen_offset_mp(var, Bbits, ctx->minfo);

        for (i = 0; i < Blen; i++)
        {
            k = Bexps[N*i + off + 0];
            check = 0;
            for (j = 1; j < wpf; j++)
                check |= Bexps[N*i + off + j];

            if (check != 0 || (slong) k < 0)
                return 0;

            fmpz_poly_set_coeff_fmpz(A, k, Bcoeffs + i);
        }
        return 1;
    }
}

void _fmpz_mpoly_set_fmpz_poly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz * Bcoeffs,
    slong Blen,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);
    slong i, Alen;
    ulong * genexp;
    TMP_INIT;

    TMP_START;

    genexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    if (Abits <= FLINT_BITS)
        mpoly_gen_monomial_sp(genexp, var, Abits, ctx->minfo);
    else
        mpoly_gen_monomial_offset_mp(genexp, var, Abits, ctx->minfo);

    Alen = 2;
    for (i = 0; i < Blen; i++)
        Alen += (Bcoeffs[i] != 0);

    fmpz_mpoly_fit_length_reset_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        if (fmpz_is_zero(Bcoeffs + i))
            continue;

        FLINT_ASSERT(Alen < A->alloc);
        fmpz_set(A->coeffs + Alen, Bcoeffs + i);
        if (Abits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps + N*Alen, genexp, N, i);
        else
            mpoly_monomial_mul_ui_mp(A->exps + N*Alen, genexp, N, i);
        Alen++;
    }
    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}

void fmpz_mpoly_set_fmpz_poly(
    fmpz_mpoly_t A,
    const fmpz_poly_t B,
    slong v,
    const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;

    if (B->length < 1)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }

    bits = mpoly_gen_pow_exp_bits_required(v, B->length - 1, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    _fmpz_mpoly_set_fmpz_poly(A, bits, B->coeffs, B->length, v, ctx);
}

/*
    construct a polynomial with at least Aminbits bits from the coefficients
    Acoeffs[0], ..., Acoeffs[deg]. *** These fmpz's are cleared *** .
*/
void _fmpz_mpoly_set_fmpz_poly_one_var(
    fmpz_mpoly_t A,
    flint_bitcnt_t Aminbits,
    fmpz * Acoeffs,     /* cleared upon return */
    slong Adeg,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, Alen;
    flint_bitcnt_t Abits;

    Abits = mpoly_gen_pow_exp_bits_required(0, Adeg, ctx->minfo);
    Abits = FLINT_MAX(Abits, Aminbits);
    Abits = mpoly_fix_bits(Abits, ctx->minfo);
    fmpz_mpoly_fit_length_reset_bits(A, Adeg + 1, Abits, ctx);

    FLINT_ASSERT(Abits <= FLINT_BITS);

    Alen = 0;
    if (ctx->minfo->ord == ORD_LEX)
    {
        FLINT_ASSERT(1 == mpoly_words_per_exp(Abits, ctx->minfo));

        for (i = Adeg; i >= 0; i--)
        {
            if (fmpz_is_zero(Acoeffs + i))
                continue;

            fmpz_swap(A->coeffs + Alen, Acoeffs + i);
            A->exps[Alen] = i;
            Alen++;
            fmpz_clear(Acoeffs + i);
        }
    }
    else if (1 == mpoly_words_per_exp(Abits, ctx->minfo))
    {
        for (i = Adeg; i >= 0; i--)
        {
            if (fmpz_is_zero(Acoeffs + i))
                continue;

            fmpz_swap(A->coeffs + Alen, Acoeffs + i);
            A->exps[Alen] = i + (i << Abits);
            Alen++;
            fmpz_clear(Acoeffs + i);
        }
    }
    else
    {
        FLINT_ASSERT(2 == mpoly_words_per_exp(Abits, ctx->minfo));

        for (i = Adeg; i >= 0; i--)
        {
            if (fmpz_is_zero(Acoeffs + i))
                continue;

            fmpz_swap(A->coeffs + Alen, Acoeffs + i);
            A->exps[2*Alen + 1] = A->exps[2*Alen + 0] = i;
            Alen++;
            fmpz_clear(Acoeffs + i);
        }
    }

    _fmpz_mpoly_set_length(A, Alen, ctx);
}

#define ALLOC_PER_VAR ((FLINT_BITS+4)/3)

char *
_fmpz_mpoly_get_str_pretty(const fmpz * coeffs, const ulong * exps, slong len,
                        const char ** x_in, slong bits, const mpoly_ctx_t mctx)
{
    char * str, ** x = (char **) x_in, *xtmp;
    slong i, j, N, bound, off;
    fmpz * exponents;
    int first;
    TMP_INIT;

    if (len == 0)
    {
        str = flint_malloc(2);
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    N = mpoly_words_per_exp(bits, mctx);

    TMP_START;

    if (x == NULL)
    {
        xtmp = (char *) TMP_ALLOC(mctx->nvars * ALLOC_PER_VAR * sizeof(char));
        x = (char **) TMP_ALLOC(mctx->nvars*sizeof(char *));
        for (i = 0; i < mctx->nvars; i++)
        {
            x[i] = xtmp + i * ALLOC_PER_VAR;
            flint_sprintf(x[i], "x%wd", i + 1);
        }
    }

    bound = 1;
    for (i = 0; i < len; i++)
        bound += fmpz_sizeinbase(coeffs + i, 10) + 1;

    exponents = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(exponents + i);
    mpoly_degrees_ffmpz((fmpz *) exponents, exps, len, bits, mctx);

    for (i = 0; i < mctx->nvars; i++)
        bound += (fmpz_sizeinbase(exponents + i, 10) + strlen(x[i]) + 3)*len;

    str = flint_malloc(bound);
    off = 0;

    for (i = 0; i < len; i++)
    {
        if (fmpz_sgn(coeffs + i) > 0 && i != 0)
            str[off++] = '+';
        if (coeffs[i] == -WORD(1))
            str[off++] = '-';
        if (coeffs[i] != WORD(1) && coeffs[i] != -WORD(1))
        {
            if (!COEFF_IS_MPZ(coeffs[i]))
                off += flint_sprintf(str + off, "%wd", coeffs[i]);
            else
                off += gmp_sprintf(str + off, "%Zd", COEFF_TO_PTR(coeffs[i]));
        }

        mpoly_get_monomial_ffmpz(exponents, exps + N*i, bits, mctx);

        first = 1;

        for (j = 0; j < mctx->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, WORD(1));
            if (cmp > 0)
            {
                if (!first || (coeffs[i] != WORD(1) && coeffs[i] != -WORD(1)))
                    off += flint_sprintf(str + off, "*");
                off += flint_sprintf(str + off, "%s^", x[j]);

                if (!COEFF_IS_MPZ(exponents[j]))
                    off += flint_sprintf(str + off, "%wd", exponents[j]);
                else
                    off += gmp_sprintf(str + off, "%Zd", COEFF_TO_PTR(exponents[j]));

                first = 0;
            }
            else if (cmp == 0)
            {
                if (!first || (coeffs[i] != WORD(1) && coeffs[i] != -WORD(1)))
                    off += flint_sprintf(str + off, "*");
                off += flint_sprintf(str + off, "%s", x[j]);
                first = 0;
            }
        }

        if (mpoly_monomial_is_zero(exps + i*N, N) && (coeffs[i] == WORD(1) || coeffs[i] == -WORD(1)))
            off += flint_sprintf(str + off, "1");
    }

    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;

    return str;
}

char *
fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx)
{
   return _fmpz_mpoly_get_str_pretty(poly->coeffs, poly->exps,
                             poly->length, x, poly->bits, ctx->minfo);
}

void fmpz_mpoly_get_term(fmpz_mpoly_t M, const fmpz_mpoly_t A,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t bits = A->bits;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_term");
    }

    fmpz_mpoly_fit_length(M, WORD(1), ctx);
    fmpz_mpoly_fit_bits(M, bits, ctx);
    M->bits = bits;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    mpoly_monomial_set(M->exps + N*0, A->exps + N*i, N);
    fmpz_set(M->coeffs + 0, A->coeffs + i);
    _fmpz_mpoly_set_length(M, 1, ctx);
}

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

void fmpz_mpoly_get_term_exp_fmpz(fmpz ** exp, const fmpz_mpoly_t A, 
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_term_exp_fmpz");
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_get_monomial_pfmpz(exp, A->exps + N*i, A->bits, ctx->minfo);
}

void fmpz_mpoly_get_term_exp_si(slong * exp, const fmpz_mpoly_t A, 
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR,
                           "Index out of range in fmpz_mpoly_get_term_exp_si");
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_get_monomial_si(exp, A->exps + N*i, A->bits, ctx->minfo);
}

void fmpz_mpoly_get_term_exp_ui(ulong * exp, const fmpz_mpoly_t A, 
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_term_exp_ui");
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_get_monomial_ui(exp, A->exps + N*i, A->bits, ctx->minfo);
}

void fmpz_mpoly_get_term_monomial(fmpz_mpoly_t M, const fmpz_mpoly_t A,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t bits = A->bits;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_term_monomial");
    }

    fmpz_mpoly_fit_length(M, WORD(1), ctx);
    fmpz_mpoly_fit_bits(M, bits, ctx);
    M->bits = bits;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    mpoly_monomial_set(M->exps + N*0, A->exps + N*i, N);
    fmpz_one(M->coeffs + 0);
    _fmpz_mpoly_set_length(M, 1, ctx);
}

slong fmpz_mpoly_get_term_var_exp_si(const fmpz_mpoly_t A, slong i,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_term_var_exp_si");
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    return mpoly_get_monomial_var_exp_si(A->exps + N*i, var, A->bits, ctx->minfo);
}

ulong fmpz_mpoly_get_term_var_exp_ui(const fmpz_mpoly_t A, slong i,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_term_var_exp_ui");
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    return mpoly_get_monomial_var_exp_ui(A->exps + N*i, var, A->bits, ctx->minfo);
}
