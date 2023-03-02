/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    assuming that the conversion can be done, set poly1 and poly1_shift
    so that poly2 is poly1 * X^poly1_shift
    TODO: handle multiprecision exponents
*/
void fmpz_mpoly_to_fmpz_poly(fmpz_poly_t poly1, slong * poly1_shift,
               const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong i, shift, off, bits, N;
    ulong k;
    slong _shift = 0, len = poly2->length;
    fmpz * coeff = poly2->coeffs;
    ulong * exp = poly2->exps;

    if (poly2->bits > FLINT_BITS)
        flint_throw(FLINT_EXPOF, "Bits too high in fmpz_mpoly_to_fmpz_poly");

    bits = poly2->bits;
    N = mpoly_words_per_exp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, var, bits, ctx->minfo);

    fmpz_poly_zero(poly1);
    if (len > 0)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        _shift = (exp[N*(len - 1)] >> shift) & mask;
        for (i = 0; i < len; i++)
        {
            k = (exp[N*i + off] >> shift) & mask;
            k -= _shift;
            FLINT_ASSERT(((slong)k) >= 0);
            fmpz_poly_set_coeff_fmpz(poly1, k, coeff + i);
        }
    }

    *poly1_shift = _shift;
}

/*
    set poly1 to poly2 * X^shift2 in the variable var
*/
void fmpz_mpoly_from_fmpz_poly(fmpz_mpoly_t poly1, const fmpz_poly_t poly2,
                           slong shift2, slong var, const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;
    slong N;
    slong k;
    slong p_len;
    fmpz * p_coeff;
    ulong * p_exp;
    slong p_alloc;
    ulong * one;
    TMP_INIT;

    TMP_START;

    bits = 1 + FLINT_BIT_COUNT(FLINT_MAX(WORD(1),
                                            shift2 + fmpz_poly_degree(poly2)));
    if (bits > FLINT_BITS)
        flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_from_fmpz_poly");
    bits = mpoly_fix_bits(bits, ctx->minfo);
    
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(one, var,bits, ctx->minfo);

    fmpz_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;

    p_coeff = poly1->coeffs;
    p_exp = poly1->exps;
    p_alloc = poly1->alloc;
    p_len = 0;
    for (k = fmpz_poly_degree(poly2); k >= 0; k--)
    {
        _fmpz_mpoly_fit_length(&p_coeff, &p_exp, &p_alloc, p_len + 1, N);
        mpoly_monomial_mul_ui(p_exp + N*p_len, one, N, k + shift2);
        fmpz_poly_get_coeff_fmpz(p_coeff + p_len, poly2, k);
        p_len += !fmpz_is_zero(p_coeff + p_len);
    }

    poly1->coeffs = p_coeff;
    poly1->exps = p_exp;
    poly1->alloc = p_alloc;
    _fmpz_mpoly_set_length(poly1, p_len, ctx);

    TMP_END;
}


/*
    set A(x_var^Bstride[var]) to B/xbar^Bshifts
    it is asserted that the conversion is correct
*/
void _fmpz_mpoly_to_fmpz_poly_deflate(
    fmpz_poly_t A,
    const fmpz_mpoly_t B,
    slong var,
    const ulong * Bshift,
    const ulong * Bstride,
    const fmpz_mpoly_ctx_t ctx)
{
    ulong mask;
    slong i, shift, off, N;
    slong len = B->length;
    fmpz * coeff = B->coeffs;
    ulong * exp = B->exps;
    ulong var_shift, var_stride;
    flint_bitcnt_t bits = B->bits;

    FLINT_ASSERT(len > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, var, bits, ctx->minfo);

    fmpz_poly_zero(A);
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    var_shift = Bshift[var];
    var_stride = Bstride[var];
    for (i = 0; i < len; i++)
    {
        ulong k = (exp[N*i + off] >> shift) & mask;
        FLINT_ASSERT(k >= var_shift);
        k -= var_shift;
        if (k != 0)
        {
            k /= var_stride;
        }
        fmpz_poly_set_coeff_fmpz(A, k, coeff + i);
    }

#if FLINT_WANT_ASSERT
    for (i = 0; i < len; i++)
    {
        slong v;
        for (v = 0; v < ctx->minfo->nvars; v++)
        {
            ulong k;
            mpoly_gen_offset_shift_sp(&off, &shift, v, bits, ctx->minfo);
            k = (exp[N*i + off] >> shift) & mask;
            FLINT_ASSERT(   (v == var && k >= Bshift[v])
                         || (v != var && k == Bshift[v]));
        }
    }
#endif
}

/*
    set A to B(x_var^Astride[var])*xbar^Ashift
    A must be packed into bits = Abits
*/
void _fmpz_mpoly_from_fmpz_poly_inflate(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_poly_t B,
    slong var,
    const ulong * Ashift,
    const ulong * Astride,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N;
    slong k;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * shiftexp;
    ulong * strideexp;
    slong Bdeg = fmpz_poly_degree(B);
    TMP_INIT;

    TMP_START;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(!fmpz_poly_is_zero(B));

    /* must have at least space for the highest exponent of var */
    FLINT_ASSERT(1 + FLINT_BIT_COUNT(Ashift[var] + Bdeg*Astride[var]) <= Abits);

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    strideexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    shiftexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_set_monomial_ui(shiftexp, Ashift, Abits, ctx->minfo);
    mpoly_gen_monomial_sp(strideexp, var, Abits, ctx->minfo);
    mpoly_monomial_mul_ui(strideexp, strideexp, N, Astride[var]);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (k = Bdeg; k >= 0; k--)
    {
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
        fmpz_poly_get_coeff_fmpz(Acoeff + Alen, B, k);
        if (!fmpz_is_zero(Acoeff + Alen))
        {
            mpoly_monomial_madd(Aexp + N*Alen, shiftexp, k, strideexp, N);
            Alen++;
        }
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}

