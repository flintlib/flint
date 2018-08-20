/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


void nmod_mpolyu_init(nmod_mpolyu_t A, mp_bitcnt_t bits,
                                                    const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void nmod_mpolyu_clear(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void nmod_mpolyu_swap(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                   const nmod_mpoly_ctx_t uctx)
{
   nmod_mpolyu_struct t = *A;
   *A = *B;
   *B = t;
}

void nmod_mpolyu_zero(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        nmod_mpoly_clear(A->coeffs + i, uctx);
        nmod_mpoly_init(A->coeffs + i, uctx);
    }
    A->length = 0;
}

int nmod_mpolyu_is_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    if (A->length != 1 || A->exps[0] != UWORD(0))
        return 0;

    return nmod_mpoly_is_one(A->coeffs + 0, uctx);
}

void nmod_mpolyu_print_pretty(const nmod_mpolyu_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        nmod_mpoly_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void nmod_mpolyu_fit_length(nmod_mpolyu_t A, slong length,
                                                   const nmod_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(nmod_mpoly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(nmod_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_mpoly_init(A->coeffs + i, uctx);
            nmod_mpoly_fit_bits(A->coeffs + i, A->bits, uctx);
            (A->coeffs + i)->bits = A->bits;
        }
        A->alloc = new_alloc;
    }
}

void nmod_mpolyu_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    nmod_mpolyu_fit_length(A, WORD(1), uctx);
    A->exps[0] = UWORD(0);
    nmod_mpoly_one(A->coeffs + 0, uctx);
    A->length = WORD(1);
}

/* if the coefficient doesn't exist, a new one is created (and set to zero) */
nmod_mpoly_struct * _nmod_mpolyu_get_coeff(nmod_mpolyu_t A,
                                        ulong pow, const nmod_mpoly_ctx_t uctx)
{
    slong i, j;
    nmod_mpoly_struct * xk;

    for (i = 0; i < A->length && A->exps[i] >= pow; i++)
    {
        if (A->exps[i] == pow) 
        {
            return A->coeffs + i;
        }
    }

    nmod_mpolyu_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        nmod_mpoly_swap(A->coeffs + j, A->coeffs + j - 1, uctx);
    }
    
    A->length++;
    A->exps[i] = pow;
    xk = A->coeffs + i;
    xk->length = 0;
    FLINT_ASSERT(xk->bits == A->bits);

    return xk;
}


void nmod_mpolyu_shift_right(nmod_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}

void nmod_mpolyu_shift_left(nmod_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        A->exps[i] += s;
    }
}

void nmod_mpolyu_scalar_mul_nmod(nmod_mpolyu_t A, mp_limb_t c,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            (A->coeffs + i)->coeffs[j] = nmod_mul((A->coeffs + i)->coeffs[j],
                                                         c, ctx->ffinfo->mod);
        }
    }
}


void nmod_mpolyu_set(nmod_mpolyu_t A, const nmod_mpolyu_t B,
                                                   const nmod_mpoly_ctx_t uctx)
{
    slong i;
    nmod_mpoly_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Alen, Blen;

    Alen = 0;
    Blen = B->length;
    nmod_mpolyu_fit_length(A, Blen, uctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpoly_set(Acoeff + Alen, Bcoeff + i, uctx);
        Aexp[Alen++] = Bexp[i];
    }
    Alen = Blen;

    /* demote remaining coefficients */
    for (i = Alen; i < A->length; i++)
    {
        nmod_mpoly_clear(Acoeff + i, uctx);
        nmod_mpoly_init(Acoeff + i, uctx);
    }
    A->length = Alen;
}



void nmod_mpoly_cvtfrom_poly_notmain(nmod_mpoly_t A, nmod_poly_t a,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    ulong * oneexp;
    slong offset, shift;
    slong N;
    TMP_INIT;
    TMP_START;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    oneexp = (ulong *)TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(oneexp, &offset, &shift, var,
                                                       N, A->bits, ctx->minfo);

    nmod_mpoly_fit_length(A, nmod_poly_length(a), ctx);

    k = 0;
    for (i = nmod_poly_length(a) - 1; i >= 0; i--)
    {
        mp_limb_t c = nmod_poly_get_coeff_ui(a, i);
        if (c != UWORD(0))
        {
            A->coeffs[k] = c;
            mpoly_monomial_mul_si(A->exps + N*k, oneexp, N, i);
            k++;
        }
    }
    A->length = k;
    TMP_END;
}
/*
    Set "A" to "a" where "a" is a polynomial in a non-main variable "var"
*/
void nmod_mpolyu_cvtfrom_poly_notmain(nmod_mpolyu_t A, nmod_poly_t a,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyu_fit_length(A, 1, ctx);
    A->exps[0] = 0;
    nmod_mpoly_cvtfrom_poly_notmain(A->coeffs + 0, a, var, ctx);
    A->length = !nmod_mpoly_is_zero(A->coeffs + 0, ctx);
}



/*
    Assuming that "A" depends only on the main variable,
    convert it to a poly "a".
*/
void nmod_mpolyu_cvtto_poly(nmod_poly_t a, nmod_mpolyu_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_poly_zero(a);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((A->coeffs + i)->length == 1);
        FLINT_ASSERT(mpoly_monomial_is_zero((A->coeffs + i)->exps, 
                      mpoly_words_per_exp((A->coeffs + i)->bits, ctx->minfo)));
        nmod_poly_set_coeff_ui(a, A->exps[i], (A->coeffs + i)->coeffs[0]);
    }
}

/*
    Convert a poly "a" to "A" in the main variable,
*/
void nmod_mpolyu_cvtfrom_poly(nmod_mpolyu_t A, nmod_poly_t a,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    nmod_mpolyu_zero(A, ctx);
    k = 0;
    for (i = nmod_poly_length(a) - 1; i >= 0; i--)
    {
        mp_limb_t c = nmod_poly_get_coeff_ui(a, i);
        if (c != UWORD(0))
        {
            nmod_mpolyu_fit_length(A, k + 1, ctx);
            A->exps[k] = i;
            nmod_mpoly_fit_length(A->coeffs + k, 1, ctx);
            nmod_mpoly_fit_bits(A->coeffs + k, A->bits, ctx);
            (A->coeffs + k)->bits = A->bits;
            (A->coeffs + k)->coeffs[0] = c;
            (A->coeffs + k)->length = 1;
            mpoly_monomial_zero((A->coeffs + k)->exps + N*0, N);
            k++;
        }
    }
    A->length = k;
}


void nmod_mpoly_to_mpolyu_perm(nmod_mpolyu_t A, const nmod_mpoly_t B,
                              const slong * perm, const nmod_mpoly_ctx_t uctx,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong n = ctx->minfo->nvars;
    slong N, NA;
    nmod_mpoly_struct * Ac;
    fmpz * uexps;
    fmpz * exps;
    TMP_INIT;

    FLINT_ASSERT(uctx->minfo->nvars == n - 1);

    TMP_START;

    uexps = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    exps  = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    for (i = 0; i < n; i++)
    {
        fmpz_init(uexps + i);
        fmpz_init(exps + i);
    }

    A->bits = B->bits;
    nmod_mpolyu_zero(A, uctx);

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    NA = mpoly_words_per_exp(A->bits, uctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ffmpz(exps, B->exps + N*j, B->bits, ctx->minfo);

        for (i = 0; i < n; i++)
        {
            fmpz_swap(uexps + i, exps + perm[i]);
        }
        Ac = _nmod_mpolyu_get_coeff(A, uexps[n - 1], uctx);

        FLINT_ASSERT(Ac->bits == B->bits);

        nmod_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        Ac->coeffs[Ac->length] = B->coeffs[j];
        mpoly_set_monomial_ffmpz(Ac->exps + NA*Ac->length, uexps, A->bits,
                                                                  uctx->minfo);
        Ac->length++;
    }

    for (i = 0; i < A->length; i++)
    {
        nmod_mpoly_sort_terms(A->coeffs + i, uctx);
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
void nmod_mpoly_from_mpolyu_perm(nmod_mpoly_t A,
                       const nmod_mpolyu_t B, int keepbits, const slong * perm,
                       const nmod_mpoly_ctx_t uctx, const nmod_mpoly_ctx_t ctx)

{
    slong i, j, k, bits;
    slong NB, N;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    slong n = ctx->minfo->nvars;
    fmpz * texps;
    fmpz * uexps;
    fmpz * exps;
    TMP_INIT;

    FLINT_ASSERT(uctx->minfo->nvars == n - 1);

    if (B->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
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
        nmod_mpoly_struct * pi = B->coeffs + i;
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

    nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        nmod_mpoly_struct * Bc = B->coeffs + i;
        _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, N);
        NB = mpoly_words_per_exp(Bc->bits, uctx->minfo);
        for (j = 0; j < Bc->length; j++)
        {
            Acoeff[Alen + j] = Bc->coeffs[j];
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
    _nmod_mpoly_set_length(A, Alen, ctx);

    for (i = 0; i < n; i++)
    {
        fmpz_clear(texps + i);
        fmpz_clear(uexps + i);
        fmpz_clear(exps + i);
    }

    nmod_mpoly_sort_terms(A, ctx);
    TMP_END;
}




void nmod_mpolyu_msub(nmod_mpolyu_t R, nmod_mpolyu_t A, nmod_mpolyu_t B,
                           nmod_mpoly_t c, slong e, const nmod_mpoly_ctx_t ctx)
{
    slong i, j, k;
    nmod_mpoly_t T;

    nmod_mpolyu_fit_length(R, A->length + B->length, ctx);

    nmod_mpoly_init(T, ctx);

    i = j = k = 0;
    while (i < A->length || j < B->length)
    {
        if (i < A->length && (j >= B->length || A->exps[i] > B->exps[j] + e))
        {
            /* only A ok */
            nmod_mpoly_set(R->coeffs + k, A->coeffs + i, ctx);
            R->exps[k] = A->exps[i];
            k++;
            i++;
        }
        else if (j < B->length && (i >= A->length || B->exps[j] + e > A->exps[i]))
        {
            /* only B ok */
            nmod_mpoly_mul(R->coeffs + k, B->coeffs + j, c, ctx);
            nmod_mpoly_neg(R->coeffs + k, R->coeffs + k, ctx);
            R->exps[k] = B->exps[j] + e;
            k++;
            j++;
        }
        else if (i < A->length && j < B->length && (A->exps[i] == B->exps[j] + e))
        {
            nmod_mpoly_mul(T, B->coeffs + j, c, ctx);
            nmod_mpoly_sub(R->coeffs + k, A->coeffs + i, T, ctx);
            R->exps[k] = A->exps[i];
            k += !nmod_mpoly_is_zero(R->coeffs + k, ctx);
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }

    nmod_mpoly_clear(T, ctx);
    R->length = k;
}

int nmod_mpolyu_divides(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int ret = 0;
    nmod_mpolyu_t P, R;
    nmod_mpoly_t t;
    nmod_mpoly_init(t, ctx);
    nmod_mpolyu_init(P, A->bits, ctx);
    nmod_mpolyu_init(R, A->bits, ctx);
    nmod_mpolyu_set(R, A, ctx);

    FLINT_ASSERT(B->length > 0);

    while (R->length > 0)
    {
        if (R->exps[0] < B->exps[0])
            goto done;

        if (!nmod_mpoly_divides(t, R->coeffs + 0, B->coeffs + 0, ctx))
            goto done;
        nmod_mpolyu_msub(P, R, B, t, R->exps[0] - B->exps[0], ctx);
        nmod_mpolyu_swap(P, R, ctx);
    }
    ret = 1;

done:
    nmod_mpoly_clear(t, ctx);
    nmod_mpolyu_clear(P, ctx);
    nmod_mpolyu_clear(R, ctx);
    return ret;
}


/*
    A = B / c and preserve the bit packing
*/
void nmod_mpolyu_divexact_mpoly(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                    nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    mp_bitcnt_t exp_bits;
    nmod_mpoly_struct * poly1, * poly2, * poly3;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == c->bits);

    nmod_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        poly1 = A->coeffs + i;
        poly2 = B->coeffs + i;
        poly3 = c;

        nmod_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
        nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        len = _nmod_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask, ctx->ffinfo);
        FLINT_ASSERT(len > 0);

        _nmod_mpoly_set_length(poly1, len, ctx);

        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
    TMP_END;
}

/*
    A = B * c and preserve the bit packing
*/
void nmod_mpolyu_mul_mpoly(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                    nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    mp_bitcnt_t exp_bits;
    nmod_mpoly_struct * poly1, * poly2, * poly3;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == c->bits);

    nmod_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        poly1 = A->coeffs + i;
        poly2 = B->coeffs + i;
        poly3 = c;

        nmod_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
        nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        len = _nmod_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps,
                       &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                        poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask, ctx->ffinfo);

        FLINT_ASSERT(len > 0);
        _nmod_mpoly_set_length(poly1, len, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
    TMP_END;
}
