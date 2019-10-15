/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"


int fq_nmod_mpolyu_is_canonical(const fq_nmod_mpolyu_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < A->length; i++)
    {
        if ((slong)(A->exps[i]) < 0)
        {
            return 0;
        }

        if (i > 0 && A->exps[i - 1] <= A->exps[i])
        {
            return 0;
        }

        if (!fq_nmod_mpoly_is_canonical(A->coeffs + i, ctx))
        {
            return 0;
        }

        if (fq_nmod_mpoly_is_zero(A->coeffs + i, ctx))
        {
            return 0;
        }
    }
    return 1;
}

void fq_nmod_mpolyu_init(fq_nmod_mpolyu_t A, flint_bitcnt_t bits,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void fq_nmod_mpolyu_clear(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fq_nmod_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

void fq_nmod_mpolyu_swap(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B)
{
   fq_nmod_mpolyu_struct t = *A;
   *A = *B;
   *B = t;
}

void fq_nmod_mpolyu_zero(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fq_nmod_mpoly_clear(A->coeffs + i, uctx);
        fq_nmod_mpoly_init(A->coeffs + i, uctx);
    }
    A->length = 0;
}

int fq_nmod_mpolyu_is_one(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
{
    if (A->length != 1 || A->exps[0] != UWORD(0))
        return 0;

    return fq_nmod_mpoly_is_one(A->coeffs + 0, uctx);
}

void fq_nmod_mpolyu_print_pretty(const fq_nmod_mpolyu_t poly,
                                const char ** x, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fq_nmod_mpoly_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fq_nmod_mpolyu_fit_length(fq_nmod_mpolyu_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fq_nmod_mpoly_struct *) flint_malloc(
                                       new_alloc*sizeof(fq_nmod_mpoly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fq_nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                       new_alloc*sizeof(fq_nmod_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fq_nmod_mpoly_init(A->coeffs + i, uctx);
            fq_nmod_mpoly_fit_bits(A->coeffs + i, A->bits, uctx);
            (A->coeffs + i)->bits = A->bits;
        }
        A->alloc = new_alloc;
    }
}

void fq_nmod_mpolyu_one(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
{
    fq_nmod_mpolyu_fit_length(A, 1, uctx);
    A->exps[0] = 0;
    fq_nmod_mpoly_one(A->coeffs + 0, uctx);
    A->length = 1;
}

/* if the coefficient doesn't exist, a new one is created (and set to zero) */
fq_nmod_mpoly_struct * _fq_nmod_mpolyu_get_coeff(fq_nmod_mpolyu_t A,
                                     ulong pow, const fq_nmod_mpoly_ctx_t uctx)
{
    slong i, j;
    fq_nmod_mpoly_struct * xk;

    for (i = 0; i < A->length && A->exps[i] >= pow; i++)
    {
        if (A->exps[i] == pow) 
        {
            return A->coeffs + i;
        }
    }

    fq_nmod_mpolyu_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        fq_nmod_mpoly_swap(A->coeffs + j, A->coeffs + j - 1, uctx);
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
void _fq_nmod_mpoly_to_mpolyu_perm_deflate(
                fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx,
                const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx,
                const slong * perm, const ulong * shift, const ulong * stride)
{
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fq_nmod_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 1 <= n);

    TMP_START;

    FLINT_ASSERT(fq_nmod_mpoly_is_canonical(B, ctx));

    uexps = (ulong *) TMP_ALLOC((m + 1)*sizeof(fmpz));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(fmpz));

    fq_nmod_mpolyu_zero(A, uctx);

    NA = mpoly_words_per_exp(A->bits, uctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        for (k = 0; k < m + 1; k++)
        {
            l = perm[k];
            FLINT_ASSERT(stride[l] != UWORD(0));
            FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == UWORD(0));
            uexps[k] = (Bexps[l] - shift[l]) / stride[l];
        }
        Ac = _fq_nmod_mpolyu_get_coeff(A, uexps[0], uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        fq_nmod_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        fq_nmod_set(Ac->coeffs + Ac->length, B->coeffs + j, ctx->fqctx);
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 1, A->bits, uctx->minfo);
        Ac->length++;
    }

    for (i = 0; i < A->length; i++)
    {
        fq_nmod_mpoly_sort_terms(A->coeffs + i, uctx);
    }

    TMP_END;

    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(A, uctx));
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
void _fq_nmod_mpoly_from_mpolyu_perm_inflate(
       fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx,
           const fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t uctx,
                const slong * perm, const ulong * shift, const ulong * stride)
{
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, k, l;
    slong NA, NB;
    slong Alen;
    fq_nmod_struct * Acoeff;
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

    fq_nmod_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_struct * Bc = B->coeffs + i;
        _fq_nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA, uctx->fqctx);
        FLINT_ASSERT(Bc->bits == B->bits);

        for (j = 0; j < Bc->length; j++)
        {
            fq_nmod_set(Acoeff + Alen + j, Bc->coeffs + j, ctx->fqctx);
            mpoly_get_monomial_ui(uexps + 1, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[0] = B->exps[i];
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k < m + 1; k++)
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
    _fq_nmod_mpoly_set_length(A, Alen, ctx);

    fq_nmod_mpoly_sort_terms(A, ctx);
    TMP_END;
}


void fq_nmod_mpolyu_shift_right(fq_nmod_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}

void fq_nmod_mpolyu_shift_left(fq_nmod_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        A->exps[i] += s;
    }
}

void fq_nmod_mpolyu_scalar_mul_fq_nmod(fq_nmod_mpolyu_t A, fq_nmod_t c,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fq_nmod_mul((A->coeffs + i)->coeffs + j,
                        (A->coeffs + i)->coeffs + j, c, ctx->fqctx);
        }
    }
}


void fq_nmod_mpolyu_set(fq_nmod_mpolyu_t A, const fq_nmod_mpolyu_t B,
                                                const fq_nmod_mpoly_ctx_t uctx)
{
    slong i;
    fq_nmod_mpoly_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Alen, Blen;

    Alen = 0;
    Blen = B->length;
    fq_nmod_mpolyu_fit_length(A, Blen, uctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        fq_nmod_mpoly_set(Acoeff + i, Bcoeff + i, uctx);
        Aexp[i] = Bexp[i];
    }
    Alen = Blen;

    /* demote remaining coefficients */
    for (i = Alen; i < A->length; i++)
    {
        fq_nmod_mpoly_clear(Acoeff + i, uctx);
        fq_nmod_mpoly_init(Acoeff + i, uctx);
    }
    A->length = Alen;
}



void fq_nmod_mpoly_cvtfrom_poly_notmain(fq_nmod_mpoly_t A, fq_nmod_poly_t a,
                                      slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    ulong * oneexp;
    slong N;
    TMP_INIT;
    TMP_START;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    oneexp = (ulong *)TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(oneexp, var, A->bits, ctx->minfo);

    fq_nmod_mpoly_fit_length(A, fq_nmod_poly_length(a, ctx->fqctx), ctx);

    k = 0;
    for (i = fq_nmod_poly_length(a, ctx->fqctx) - 1; i >= 0; i--)
    {
        mpoly_monomial_mul_ui(A->exps + N*k, oneexp, N, i);
        fq_nmod_poly_get_coeff(A->coeffs + k, a, i, ctx->fqctx);
        k += !fq_nmod_is_zero(A->coeffs + k, ctx->fqctx);
    }
    A->length = k;
    TMP_END;
}
/*
    Set "A" to "a" where "a" is a polynomial in a non-main variable "var"
*/
void fq_nmod_mpolyu_cvtfrom_poly_notmain(fq_nmod_mpolyu_t A, fq_nmod_poly_t a,
                                      slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpolyu_fit_length(A, 1, ctx);
    A->exps[0] = 0;
    fq_nmod_mpoly_cvtfrom_poly_notmain(A->coeffs + 0, a, var, ctx);
    A->length = !fq_nmod_mpoly_is_zero(A->coeffs + 0, ctx);
}



/*
    Assuming that "A" depends only on the main variable,
    convert it to a poly "a".
*/
void fq_nmod_mpolyu_cvtto_poly(fq_nmod_poly_t a, fq_nmod_mpolyu_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(A, ctx));

    fq_nmod_poly_zero(a, ctx->fqctx);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((A->coeffs + i)->length == 1);
        FLINT_ASSERT(mpoly_monomial_is_zero((A->coeffs + i)->exps, 
                      mpoly_words_per_exp((A->coeffs + i)->bits, ctx->minfo)));
        fq_nmod_poly_set_coeff(a, A->exps[i], (A->coeffs + i)->coeffs + 0,
                                                                   ctx->fqctx);
    }
}

/*
    Convert a poly "a" to "A" in the main variable,
*/
void fq_nmod_mpolyu_cvtfrom_poly(fq_nmod_mpolyu_t A, fq_nmod_poly_t a,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    fq_nmod_t c;
    fq_nmod_init(c, ctx->fqctx);

    fq_nmod_mpolyu_zero(A, ctx);
    k = 0;
    for (i = fq_nmod_poly_length(a, ctx->fqctx) - 1; i >= 0; i--)
    {
        fq_nmod_poly_get_coeff(c, a, i, ctx->fqctx);
        if (!fq_nmod_is_zero(c, ctx->fqctx))
        {
            fq_nmod_mpolyu_fit_length(A, k + 1, ctx);
            A->exps[k] = i;
            fq_nmod_mpoly_fit_length(A->coeffs + k, 1, ctx);
            fq_nmod_mpoly_fit_bits(A->coeffs + k, A->bits, ctx);
            (A->coeffs + k)->bits = A->bits;
            fq_nmod_set((A->coeffs + k)->coeffs + 0, c, ctx->fqctx);
            (A->coeffs + k)->length = 1;
            mpoly_monomial_zero((A->coeffs + k)->exps + N*0, N);
            k++;
        }
    }
    A->length = k;

    fq_nmod_clear(c, ctx->fqctx);
}


void fq_nmod_mpolyu_msub(fq_nmod_mpolyu_t R, fq_nmod_mpolyu_t A,
                             fq_nmod_mpolyu_t B, fq_nmod_mpoly_t c, slong e,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j, k;
    fq_nmod_mpoly_t T;

    fq_nmod_mpolyu_fit_length(R, A->length + B->length, ctx);

    fq_nmod_mpoly_init(T, ctx);

    i = j = k = 0;
    while (i < A->length || j < B->length)
    {
        if (i < A->length && (j >= B->length || A->exps[i] > B->exps[j] + e))
        {
            /* only A ok */
            fq_nmod_mpoly_set(R->coeffs + k, A->coeffs + i, ctx);
            R->exps[k] = A->exps[i];
            k++;
            i++;
        }
        else if (j < B->length && (i >= A->length || B->exps[j] + e > A->exps[i]))
        {
            /* only B ok */
            fq_nmod_mpoly_mul_johnson(R->coeffs + k, B->coeffs + j, c, ctx);
            fq_nmod_mpoly_neg(R->coeffs + k, R->coeffs + k, ctx);
            R->exps[k] = B->exps[j] + e;
            k++;
            j++;
        }
        else if (i < A->length && j < B->length && (A->exps[i] == B->exps[j] + e))
        {
            fq_nmod_mpoly_mul_johnson(T, B->coeffs + j, c, ctx);
            fq_nmod_mpoly_sub(R->coeffs + k, A->coeffs + i, T, ctx);
            R->exps[k] = A->exps[i];
            k += !fq_nmod_mpoly_is_zero(R->coeffs + k, ctx);
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }

    fq_nmod_mpoly_clear(T, ctx);
    R->length = k;
}

int fq_nmod_mpolyu_divides(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    int ret = 0;
    fq_nmod_mpolyu_t P, R;
    fq_nmod_mpoly_t t;
    fq_nmod_mpoly_init(t, ctx);
    fq_nmod_mpolyu_init(P, A->bits, ctx);
    fq_nmod_mpolyu_init(R, A->bits, ctx);
    fq_nmod_mpolyu_set(R, A, ctx);

    FLINT_ASSERT(B->length > 0);

    while (R->length > 0)
    {
        if (R->exps[0] < B->exps[0])
            goto done;

        if (!fq_nmod_mpoly_divides_monagan_pearce(t, R->coeffs + 0, B->coeffs + 0, ctx))
            goto done;

        fq_nmod_mpolyu_msub(P, R, B, t, R->exps[0] - B->exps[0], ctx);
        fq_nmod_mpolyu_swap(P, R);
    }
    ret = 1;

done:
    fq_nmod_mpoly_clear(t, ctx);
    fq_nmod_mpolyu_clear(P, ctx);
    fq_nmod_mpolyu_clear(R, ctx);
    return ret;
}


/*
    A = B / c and preserve the bit packing
*/
void fq_nmod_mpolyu_divexact_mpoly(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                              fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    flint_bitcnt_t exp_bits;
    fq_nmod_mpoly_struct * poly1, * poly2, * poly3;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == c->bits);

    fq_nmod_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        poly1 = A->coeffs + i;
        poly2 = B->coeffs + i;
        poly3 = c;

        fq_nmod_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
        fq_nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        len = _fq_nmod_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask, ctx->fqctx);
        FLINT_ASSERT(len > 0);

        _fq_nmod_mpoly_set_length(poly1, len, ctx);

        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
    TMP_END;
}

/*
    A = B * c and preserve the bit packing
*/
void fq_nmod_mpolyu_mul_mpoly(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                              fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    flint_bitcnt_t exp_bits;
    fq_nmod_mpoly_struct * poly1, * poly2, * poly3;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == c->bits);

    fq_nmod_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        poly1 = A->coeffs + i;
        poly2 = B->coeffs + i;
        poly3 = c;

        fq_nmod_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
        fq_nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        len = _fq_nmod_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps,
                      &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                       poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask, ctx->fqctx);

        FLINT_ASSERT(len > 0);
        _fq_nmod_mpoly_set_length(poly1, len, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
    TMP_END;
}

int fq_nmod_mpolyu_content_mpoly(
    fq_nmod_mpoly_t g,
    const fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    int success;
    flint_bitcnt_t bits = A->bits;

    FLINT_ASSERT(g->bits == bits);

    if (A->length < 2)
    {
        if (A->length == 0)
        {
            fq_nmod_mpoly_zero(g, ctx);
        }
        else
        {
            fq_nmod_mpoly_make_monic(g, A->coeffs + 0, ctx);
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
        j = 1;

    success = _fq_nmod_mpoly_gcd(g, bits, A->coeffs + 0, A->coeffs + j, ctx);
    if (!success)
        return 0;

    for (i = 1; i < A->length; i++)
    {
        if (i == j)
            continue;

        success = _fq_nmod_mpoly_gcd(g, bits, g, A->coeffs + i, ctx);
        FLINT_ASSERT(g->bits == bits);
        if (!success)
            return 0;
    }

    return 1;
}
