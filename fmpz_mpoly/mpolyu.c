/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"


int fmpz_mpolyu_is_canonical(const fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (A->length > A->alloc)
    {
        return 0;
    }

    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_mpoly_is_canonical(A->coeffs + i, ctx))
        {
            return 0;
        }

        if (i > 0 && A->exps[i - 1] <= A->exps[i])
        {
            return 0;
        }
    }

    return 1;
}

void fmpz_mpolyu_init(fmpz_mpolyu_t A, mp_bitcnt_t bits,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void fmpz_mpolyu_clear(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void fmpz_mpolyu_swap(fmpz_mpolyu_t A, fmpz_mpolyu_t B,
                                                   const fmpz_mpoly_ctx_t uctx)
{
   fmpz_mpolyu_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mpolyu_zero(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fmpz_mpoly_clear(A->coeffs + i, uctx);
        fmpz_mpoly_init(A->coeffs + i, uctx);
    }
    A->length = 0;
}



void fmpz_mpolyu_print_pretty(
    const fmpz_mpolyu_t poly,
    const char ** x,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fmpz_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")*X^%wd", poly->exps[i]);
    }
}

void fmpz_mpolyuu_print_pretty(
    const fmpz_mpolyu_t poly,
    const char ** x,
    slong nmainvars,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/nmainvars);

    if (poly->length == 0)
        flint_printf("0");

    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fmpz_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")");
        for (j = nmainvars - 1; j >= 0; j--)
        {
            flint_printf("*X%wd^%wd", nmainvars - 1 - j,
                    mask & (poly->exps[i] >> (FLINT_BITS/nmainvars*j)));
        }
    }
}


void fmpz_mpolyu_fit_length(fmpz_mpolyu_t A, slong length,
                                                   const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpoly_init(A->coeffs + i, uctx);
            fmpz_mpoly_fit_bits(A->coeffs + i, A->bits, uctx);
            (A->coeffs + i)->bits = A->bits;
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mpolyu_one(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx)
{
    fmpz_mpolyu_fit_length(A, WORD(1), uctx);
    A->exps[0] = UWORD(0);
    fmpz_mpoly_one(A->coeffs + 0, uctx);
    A->length = WORD(1);
}


void fmpz_mpolyu_set(fmpz_mpolyu_t A, const fmpz_mpolyu_t B,
                                                   const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    fmpz_mpoly_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Alen, Blen;

    Alen = 0;
    Blen = B->length;
    fmpz_mpolyu_fit_length(A, Blen, uctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        fmpz_mpoly_set(Acoeff + Alen, Bcoeff + i, uctx);
        Aexp[Alen++] = Bexp[i];
    }
    Alen = Blen;

    /* demote remaining coefficients */
    for (i = Alen; i < A->length; i++)
    {
        fmpz_mpoly_clear(Acoeff + i, uctx);
        fmpz_mpoly_init(Acoeff + i, uctx);
    }
    A->length = Alen;
}


/* if the coefficient doesn't exist, a new one is created (and set to zero) */
fmpz_mpoly_struct * _fmpz_mpolyu_get_coeff(fmpz_mpolyu_t A,
                             ulong pow, const fmpz_mpoly_ctx_t uctx)
{
    slong i, j;
    fmpz_mpoly_struct * xk;

    for (i = 0; i < A->length && A->exps[i] >= pow; i++)
    {
        if (A->exps[i] == pow) 
        {
            return A->coeffs + i;
        }
    }

    fmpz_mpolyu_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        fmpz_mpoly_swap(A->coeffs + j, A->coeffs + j - 1, uctx);
    }

    A->length++;
    A->exps[i] = pow;
    xk = A->coeffs + i;
    xk->length = 0;
    FLINT_ASSERT(xk->bits == A->bits);

    return xk;
}


/*
    Convert B to A using the variable permutation perm.
    The uctx should be the context of the coefficients of A.
    The ctx should be the context of B.

    operation on each term:

    for 0 <= k <= m
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]
*/
void fmpz_mpoly_to_mpolyu_perm_deflate(fmpz_mpolyu_t A, const fmpz_mpoly_t B,
               const slong * perm, const ulong * shift, const ulong * stride,
                       const fmpz_mpoly_ctx_t uctx, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fmpz_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 1 <= n);

    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 1)*sizeof(fmpz));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(fmpz));

    fmpz_mpolyu_zero(A, uctx);

    NA = mpoly_words_per_exp(A->bits, uctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        for (k = 0; k <= m; k++)
        {
            l = perm[k];
            FLINT_ASSERT(stride[l] != UWORD(0));
            FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == UWORD(0));
            uexps[k] = (Bexps[l] - shift[l]) / stride[l];
        }
        Ac = _fmpz_mpolyu_get_coeff(A, uexps[m], uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        fmpz_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps, A->bits, uctx->minfo);
        Ac->length++;
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_sort_terms(A->coeffs + i, uctx);
    }

    TMP_END;
}


/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k <= m
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void fmpz_mpoly_from_mpolyu_perm_inflate(fmpz_mpoly_t A, mp_bitcnt_t Abits,
                                                        const fmpz_mpolyu_t B,
                const slong * perm, const ulong * shift, const ulong * stride,
                       const fmpz_mpoly_ctx_t uctx, const fmpz_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, k, l;
    slong NA, NB;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * uexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 1 <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 1)*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, uctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA);
        FLINT_ASSERT(Bc->bits == B->bits);

        for (j = 0; j < Bc->length; j++)
        {
            fmpz_set(Acoeff + Alen + j, Bc->coeffs + j);
            mpoly_get_monomial_ui(uexps, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[m] = B->exps[i];
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k <= m; k++)
            {
                l = perm[k];
                Aexps[l] += stride[l]*uexps[k];
            }
            mpoly_set_monomial_ui(Aexp + NA*(Alen + j), Aexps, Abits, ctx->minfo);
        }
        Alen += Bc->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    fmpz_mpoly_sort_terms(A, ctx);
    TMP_END;
}


/*
    Convert B to A using the variable permutation perm.
    The uctx should be the context of the coefficients of A.
    The ctx should be the context of B.

    operation on each term:

    for 0 <= k < m + 2
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]

    the most significant main variable uses Aexp[0]
    the least significant main variable uses Aexp[1]
    the coefficients of A use variables Aexp[2], ..., Aexp[m + 1]
*/
void fmpz_mpoly_to_mpolyuu_perm_deflate(
    fmpz_mpolyu_t A,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const fmpz_mpoly_ctx_t uctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fmpz_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);

    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    fmpz_mpolyu_zero(A, uctx);

    NA = mpoly_words_per_exp(A->bits, uctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        for (k = 0; k < m + 2; k++)
        {
            l = perm[k];
            FLINT_ASSERT(stride[l] != UWORD(0));
            FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == UWORD(0));
            uexps[k] = (Bexps[l] - shift[l]) / stride[l];
        }
        FLINT_ASSERT(FLINT_BIT_COUNT(uexps[0]) < FLINT_BITS/2);
        FLINT_ASSERT(FLINT_BIT_COUNT(uexps[1]) < FLINT_BITS/2);
        Ac = _fmpz_mpolyu_get_coeff(A, (uexps[0] << (FLINT_BITS/2)) + uexps[1], uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        fmpz_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 2, A->bits, uctx->minfo);
        Ac->length++;
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_sort_terms(A->coeffs + i, uctx);
    }

    TMP_END;
}


/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k < m + 2
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void fmpz_mpoly_from_mpolyuu_perm_inflate( /* only for 2 main vars */
    fmpz_mpoly_t A,
    mp_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mpolyu_t B,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const fmpz_mpoly_ctx_t uctx)
{
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, k, l;
    slong NA, NB;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * uexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, uctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA);
        FLINT_ASSERT(Bc->bits == B->bits);

        for (j = 0; j < Bc->length; j++)
        {
            fmpz_set(Acoeff + Alen + j, Bc->coeffs + j);
            mpoly_get_monomial_ui(uexps + 2, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[0] = B->exps[i] >> (FLINT_BITS/2);
            uexps[1] = B->exps[i] & ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2));
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k < m + 2; k++)
            {
                l = perm[k];
                Aexps[l] += stride[l]*uexps[k];
            }
            mpoly_set_monomial_ui(Aexp + NA*(Alen + j), Aexps, Abits, ctx->minfo);
        }
        Alen += Bc->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    fmpz_mpoly_sort_terms(A, ctx);
    TMP_END;
}



/*
    Convert B to A using the variable permutation perm.
    The uctx should be the context of the coefficients of A.
    The ctx should be the context of B.
    The uctx should have one fewer variable than ctx.
*/
void fmpz_mpoly_to_mpolyu_perm(fmpz_mpolyu_t A, const fmpz_mpoly_t B,
                              const slong * perm, const fmpz_mpoly_ctx_t uctx,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong n = ctx->minfo->nvars;
    slong N, NA;
    fmpz * uexps;
    fmpz * exps;
    fmpz_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(uctx->minfo->nvars == n - 1);

    TMP_START;

    uexps = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    exps  = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    for (i = 0; i < n; i++)
    {
        fmpz_init(uexps + i);
        fmpz_init(exps + i);
    }

    fmpz_mpolyu_zero(A, uctx);

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    NA = mpoly_words_per_exp(A->bits, uctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ffmpz(exps, B->exps + N*j, B->bits, ctx->minfo);

        for (i = 0; i < n; i++)
        {
            fmpz_swap(uexps + i, exps + perm[i]);
        }
        Ac = _fmpz_mpolyu_get_coeff(A, uexps[n - 1], uctx);

        FLINT_ASSERT(Ac->bits == A->bits);

        fmpz_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
        mpoly_set_monomial_ffmpz(Ac->exps + NA*Ac->length, uexps, A->bits, uctx->minfo);
        Ac->length++;
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_sort_terms(A->coeffs + i, uctx);
    }

    for (i = 0; i < n; i++)
    {
        fmpz_clear(uexps + i);
        fmpz_clear(exps + i);
    }

    TMP_END;
}



/*
    Convert B to A using the variable permutation vector perm.
    If keepbits != 0, A is constructed with the same bit count as B.
    If keepbits == 0, A is constructed with the default bit count.
*/
void fmpz_mpoly_from_mpolyu_perm(fmpz_mpoly_t A,
                       const fmpz_mpolyu_t B, int keepbits, const slong * perm,
                       const fmpz_mpoly_ctx_t uctx, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k;
    slong N, NB;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    slong n = ctx->minfo->nvars;
    fmpz * texps;
    fmpz * uexps;
    fmpz * exps;
    mp_bitcnt_t bits;
    TMP_INIT;

    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(uctx->minfo->nvars == n - 1);

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;        
    }

    TMP_START;

    /* find bits required to represent result */
    texps = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    uexps = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    exps  = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    for (i = 0; i < n; i++)
    {
        fmpz_init(texps + i);
        fmpz_init(uexps + i);
        fmpz_init(exps + i);
    }
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * pi = B->coeffs + i;
        mpoly_degrees_ffmpz(texps, pi->exps, pi->length, pi->bits, uctx->minfo);
        _fmpz_vec_max_inplace(uexps, texps, n - 1);
    }
    fmpz_set_ui(uexps + n - 1, B->exps[0]);
    for (i = 0; i < n; i++)
    {
        fmpz_swap(uexps + i, exps + perm[i]);
    }
    bits = mpoly_exp_bits_required_ffmpz(exps, ctx->minfo);

    if (keepbits)
    {
        /* the bit count of B should be sufficient to represent A */
        FLINT_ASSERT(bits <= B->bits);
        bits = B->bits;
    } else {
        /* upgrade the bit count to the default */
        bits = FLINT_MAX(MPOLY_MIN_BITS, bits);
        bits = mpoly_fix_bits(bits, ctx->minfo);
    }

    N = mpoly_words_per_exp(bits, ctx->minfo);

    fmpz_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, N);
        NB = mpoly_words_per_exp(Bc->bits, uctx->minfo);
        for (j = 0; j < Bc->length; j++)
        {
            fmpz_set(Acoeff + Alen + j, Bc->coeffs + j);
            mpoly_get_monomial_ffmpz(uexps, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            fmpz_set_ui(uexps + n - 1, B->exps[i]);

            for (k = 0; k < n; k++)
            {
                fmpz_swap(uexps + k, exps + perm[k]);
            }

            mpoly_set_monomial_ffmpz(Aexp + N*(Alen + j), exps, bits, ctx->minfo);
        }
        Alen += (B->coeffs + i)->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    for (i = 0; i < n; i++)
    {
        fmpz_clear(texps + i);
        fmpz_clear(uexps + i);
        fmpz_clear(exps + i);
    }

    fmpz_mpoly_sort_terms(A, ctx);
    TMP_END;
}


void fmpz_mpolyu_fmpz_content(fmpz_t c, fmpz_mpolyu_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    fmpz_zero(c);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fmpz_gcd(c, c, (A->coeffs + i)->coeffs + j);
            if (fmpz_is_one(c))
                break;
        }
    }
}



void fmpz_mpolyu_divexact_fmpz(
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    fmpz_t c,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(!fmpz_is_zero(c));
    FLINT_ASSERT(A->bits == B->bits);
    fmpz_mpolyu_fit_length(A, B->length, ctx);

    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fmpz_mpoly_scalar_divexact_fmpz(A->coeffs + i, B->coeffs + i, c, ctx);
        FLINT_ASSERT((A->coeffs + i)->bits == B->bits);
    }
    A->length = B->length;
}


/* The bit counts of A, B and c must all agree for this division */
void fmpz_mpolyu_divexact_mpoly(
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    fmpz_mpoly_t c,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    mp_bitcnt_t exp_bits;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == c->bits);

    fmpz_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * poly1 = A->coeffs + i;
        fmpz_mpoly_struct * poly2 = B->coeffs + i;
        fmpz_mpoly_struct * poly3 = c;
        A->exps[i] = B->exps[i];

        fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
        fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        FLINT_ASSERT(poly2->length > 0);

        len = _fmpz_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask);
        FLINT_ASSERT(len > 0);
        _fmpz_mpoly_set_length(poly1, len, ctx);
    }
    A->length = B->length;

    TMP_END;
}



/*
    The bit counts of A, B and c must all agree for this multiplication
*/
void fmpz_mpolyu_mul_mpoly(fmpz_mpolyu_t A, fmpz_mpolyu_t B,
                                    fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    mp_bitcnt_t exp_bits;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == c->bits);

    fmpz_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * poly1 = A->coeffs + i;
        fmpz_mpoly_struct * poly2 = B->coeffs + i;
        fmpz_mpoly_struct * poly3 = c;
        A->exps[i] = B->exps[i];

        fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
        fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        len = _fmpz_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask);

        _fmpz_mpoly_set_length(poly1, len, ctx);

    }
    A->length = B->length;

    TMP_END;
}


void fmpz_mpolyu_shift_right(fmpz_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}


void fmpz_mpolyu_shift_left(fmpz_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((slong)(A->exps[i] + s) >= 0);
        A->exps[i] += s;
    }
}


void fmpz_mpolyu_content_fmpz(
    fmpz_t g,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    fmpz_zero(g);
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_gcd(g, g, Ac->coeffs + j);
            if (fmpz_is_one(g))
                return;
        }
    }
}


int fmpz_mpolyu_content_mpoly(
    fmpz_mpoly_t g,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    int success;
    mp_bitcnt_t bits = A->bits;

    FLINT_ASSERT(g->bits == bits);

    if (A->length < 2)
    {
        if (A->length == 0)
        {
            fmpz_mpoly_zero(g, ctx);
        }
        else
        {
            fmpz_mpoly_set(g, A->coeffs + 0, ctx);
        }

        FLINT_ASSERT(g->bits == bits);
        return 1;
    }

    j = 0;
    for (i = 1; i < A->length; i++)
    {
        if ((A->coeffs + i)->length < (A->coeffs + j)->length)
        {
            j = i;
        }
    }

    if (j == 0)
    {
        j = 1;
    }
    success = _fmpz_mpoly_gcd(g, bits, A->coeffs + 0, A->coeffs + j, ctx);
    if (!success)
    {
        return 0;
    }

    for (i = 1; i < A->length; i++)
    {
        if (i == j)
        {
            continue;
        }
        success = _fmpz_mpoly_gcd(g, bits, g, A->coeffs + i, ctx);
        FLINT_ASSERT(g->bits == bits);
        if (!success)
        {
            return 0;
        }
    }

    return 1;
}
