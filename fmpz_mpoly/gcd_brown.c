/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "nmod_mpoly.h"


void fmpz_mpolyd_swap(fmpz_mpolyd_t A, fmpz_mpolyd_t B)
{
    fmpz_mpolyd_struct t = *A;
    *A = *B;
    *B = t;
}

void fmpz_mpolyd_ctx_init(fmpz_mpolyd_ctx_t dctx, slong nvars)
{
    slong i;

    dctx->nvars = nvars;
    dctx->perm = (slong *) flint_malloc(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }
}


int fmpz_mpolyd_ctx_init_version1(fmpz_mpolyd_ctx_t dctx,
                            const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success = 0;
    slong i, j, degb_prod;
    slong * Aexps, * Bexps, * deg_bounds;
    slong nvars = ctx->minfo->nvars;
    slong * perm = dctx->perm;
    TMP_INIT;

    TMP_START;
    Aexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Bexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
        goto cleanup;
    fmpz_mpoly_degrees_si(Aexps, A, ctx);
    fmpz_mpoly_degrees_si(Bexps, B, ctx);

    deg_bounds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }

    degb_prod = 1;
    for (i = 0; i < nvars; i++)
    {
        ulong hi;
        deg_bounds[i] = FLINT_MAX(Aexps[i] + 1, Bexps[i] + 1);
        umul_ppmm(hi, degb_prod, degb_prod, deg_bounds[i]);
        if (hi != WORD(0) || degb_prod < 0)
            goto cleanup;
    }

    success = 1;
    for (i = 1; i < nvars; i++)
    {
        for (j = i; (j > 0) && (deg_bounds[j-1] < deg_bounds[j-0]); j--)
        {
            slong t1, t2;
            t1 = deg_bounds[j-1];
            t2 = deg_bounds[j-0];
            deg_bounds[j-0] = t1;
            deg_bounds[j-1] = t2;
            t1 = perm[j-1];
            t2 = perm[j-0];
            perm[j-0] = t1;
            perm[j-1] = t2;
        }
    }

cleanup:
    TMP_END;
    return success;
}


void fmpz_mpolyd_ctx_clear(fmpz_mpolyd_ctx_t dctx)
{
    flint_free(dctx->perm);
}


void fmpz_mpolyd_init(fmpz_mpolyd_t poly, slong nvars)
{
    slong i;

    poly->nvars = nvars;
    poly->degb_alloc = nvars;
    poly->deg_bounds = (slong *) flint_malloc(poly->degb_alloc*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeff_alloc = WORD(16);
    poly->coeffs = (fmpz *) flint_malloc(poly->coeff_alloc*sizeof(fmpz));
    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fmpz_init(poly->coeffs + i);
    }
}


void fmpz_mpolyd_fit_length(fmpz_mpolyd_t poly, slong len)
{
    if (poly->coeff_alloc < len) {
        slong i;
        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, len*sizeof(fmpz));
        for (i = poly->coeff_alloc; i < len; i++)
        {
            fmpz_init(poly->coeffs + i);
        }
        poly->coeff_alloc = len;
    }
}


void fmpz_mpolyd_set_nvars(fmpz_mpolyd_t poly, slong nvars)
{
    poly->nvars = nvars;
    if (poly->degb_alloc < nvars) {
        poly->deg_bounds = (slong *) flint_realloc(poly->deg_bounds, nvars*sizeof(slong));
        poly->degb_alloc = nvars;
    }
}


void fmpz_mpolyd_zero(fmpz_mpolyd_t poly)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeffs[0] = UWORD(0);
}

void fmpz_mpolyd_set_fmpz(fmpz_mpolyd_t poly, fmpz_t num)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    fmpz_set(poly->coeffs + 0, num);
}


void fmpz_mpolyd_clear(fmpz_mpolyd_t poly)
{
    slong i;

    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fmpz_clear(poly->coeffs + i);
    }

    flint_free(poly->deg_bounds);
    flint_free(poly->coeffs);
    poly->deg_bounds = NULL;
    poly->coeffs = NULL;
}


void fmpz_mpoly_convert_to_fmpz_mpolyd(
                            fmpz_mpolyd_t poly1, const fmpz_mpolyd_ctx_t dctx,
                          const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong degb_prod;
    slong i, j, N;
    slong * exps;
    const slong * perm = dctx->perm;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    fmpz_mpolyd_set_nvars(poly1, ctx->minfo->nvars);

    FLINT_ASSERT(poly2->bits <= FLINT_BITS);

    if (poly2->length == 0)
    {
        fmpz_mpolyd_zero(poly1);
        return;
    }

    TMP_START;
    exps = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));

    fmpz_mpoly_degrees_si(exps, poly2, ctx);
    degb_prod = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        poly1->deg_bounds[i] = exps[perm[i]] + 1;
        degb_prod *= poly1->deg_bounds[i];
    }

    fmpz_mpolyd_fit_length(poly1, degb_prod);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_zero(poly1->coeffs + i);
    }

    N = mpoly_words_per_exp(poly2->bits, ctx->minfo);
    for (i = 0; i < poly2->length; i++)
    {
        slong off = 0;

        mpoly_get_monomial_ui((ulong *)exps, poly2->exps + N*i, poly2->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            off = exps[perm[j]] + poly1->deg_bounds[j]*off;
        }

        fmpz_set(poly1->coeffs + off, poly2->coeffs + i);
    }

    TMP_END;
}


/*
    m is the number of variables in A
*/
void fmpz_mpoly_to_fmpz_mpolyd_perm_deflate(fmpz_mpolyd_t A, slong m,
              const fmpz_mpoly_t B, const slong * perm, const ulong * shift,
        const ulong * stride, const ulong * degree, const fmpz_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong degb_prod;
    slong i, k, l, N;
    ulong * Bexp;
    TMP_INIT;

    FLINT_ASSERT(m <= n);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(B->length > 0);

    fmpz_mpolyd_set_nvars(A, m);

    TMP_START;
    Bexp = (ulong *) TMP_ALLOC(n*sizeof(slong));

    degb_prod = WORD(1);
    for (k = 0; k < m; k++)
    {
        l = perm[k];
        FLINT_ASSERT(stride[l] != UWORD(0));
        FLINT_ASSERT((degree[l] - shift[l]) % stride[l] == UWORD(0));
        A->deg_bounds[k] = (degree[l] - shift[l])/stride[l] + 1;
        degb_prod *= A->deg_bounds[k];
        /* we should not be converting something whose dense size overflows */
        FLINT_ASSERT(degb_prod > 0);
        FLINT_ASSERT(degb_prod >= A->deg_bounds[k]);
    }

    fmpz_mpolyd_fit_length(A, degb_prod);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_zero(A->coeffs + i);
    }

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        slong off = 0;
        mpoly_get_monomial_ui(Bexp, B->exps + N*i, B->bits, ctx->minfo);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            FLINT_ASSERT(stride[l] != UWORD(0));
            FLINT_ASSERT(((Bexp[l] - shift[l]) % stride[l]) == UWORD(0));
            FLINT_ASSERT((Bexp[l] - shift[l])/stride[l] < A->deg_bounds[k]);
            off = (Bexp[l] - shift[l])/stride[l] + A->deg_bounds[k]*off;
        }
        fmpz_set(A->coeffs + off, B->coeffs + i);
    }

    TMP_END;
}


void fmpz_mpoly_from_fmpz_mpolyd_perm_inflate(fmpz_mpoly_t A,
         mp_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx, const fmpz_mpolyd_t B,
                 const slong * perm, const ulong * shift, const ulong * stride)
{
    slong off;
    slong n = ctx->minfo->nvars;
    slong m = B->nvars;
    slong Alen;
    slong i, j, l, k, N;
    slong perm_nontrivial;
    ulong topmask;
    ulong * exps, * pcurexp, * pexps;
    TMP_INIT;

    FLINT_ASSERT(m <= n);
    FLINT_ASSERT(Abits <= FLINT_BITS);

    perm_nontrivial = n - m;

    /* we are going to push back terms manually */
    Alen = 0;
    fmpz_mpoly_zero(A, ctx);
    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    N = mpoly_words_per_exp(Abits, ctx->minfo);

    TMP_START;

    /* find exponent vector for all variables in B */
    pexps = (ulong *) TMP_ALLOC(N*m*sizeof(ulong));
    for (k = 0; k < m; k++)
    {
        l = perm[k];
        perm_nontrivial |= l - k;
        mpoly_gen_monomial_sp(pexps + k*N, l, Abits, ctx->minfo);
        mpoly_monomial_mul_ui(pexps + k*N, pexps + k*N, N, stride[l]);
    }

    /* get most significant exponent in pcurexp and its vector in exps */
    pcurexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    off = WORD(1);
    for (j = 0; j < m; j++)
    {
        off *= B->deg_bounds[j];
    }
    FLINT_ASSERT(off <= B->coeff_alloc);
    off--;
    mpoly_set_monomial_ui(pcurexp, shift, Abits, ctx->minfo);
    i = off;
    for (k = m - 1; k >= 0; k--) 
    {
        exps[k] = i % B->deg_bounds[k];
        i = i / B->deg_bounds[k];
        mpoly_monomial_madd(pcurexp, pcurexp, exps[k], pexps + N*k, N);
    }

    /* scan down through the exponents */
    topmask = 0;

    for (; off >= 0; off--)
    {
        if (!fmpz_is_zero(B->coeffs + off))
        {
            _fmpz_mpoly_fit_length(&A->coeffs, &A->exps, &A->alloc, Alen + 1, N);
            fmpz_set(A->coeffs + Alen, B->coeffs + off);
            mpoly_monomial_set(A->exps + N*Alen, pcurexp, N);
            topmask |= pcurexp[N - 1];
            Alen++;
        }

        k = m - 1;
        do {
            --exps[k];
            if ((slong)(exps[k]) < WORD(0))
            {
                FLINT_ASSERT(off == 0 || k > 0);
                FLINT_ASSERT(exps[k] == -UWORD(1));
                exps[k] = B->deg_bounds[k] - 1;
                mpoly_monomial_madd(pcurexp, pcurexp, exps[k], pexps + N*k, N);
            }
            else
            {
                mpoly_monomial_sub(pcurexp, pcurexp, pexps + N*k, N);
                break;
            }
        } while (--k >= 0);
    }
    _fmpz_mpoly_set_length(A, Alen, ctx);


    /* sort the terms if needed */
    if (ctx->minfo->ord != ORD_LEX || perm_nontrivial != WORD(0))
    {
        slong msb;
        mpoly_get_cmpmask(pcurexp, N, Abits, ctx->minfo);
        if (topmask != WORD(0))
        {
            count_leading_zeros(msb, topmask);
            msb = (FLINT_BITS - 1)^msb;
        }
        else
        {
            msb = -WORD(1);
        }
        if (N == 1)
        {
            if (msb >= WORD(0))
            {
                _fmpz_mpoly_radix_sort1(A, 0, A->length,
                                                   msb, pcurexp[0], topmask);
            }
        }
        else
        {
            _fmpz_mpoly_radix_sort(A, 0, A->length,
                                        (N - 1)*FLINT_BITS + msb, N, pcurexp);
        }
    }

    TMP_END;
}


void fmpz_mpolyd_print(fmpz_mpolyd_t poly, const char ** vars,
                                                  const fmpz_mpolyd_ctx_t dctx)
{
    int first = 0;
    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < poly->nvars; j++) {
        degb_prod *= poly->deg_bounds[j];
    }

    first = 1;
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (fmpz_is_zero(poly->coeffs + i))
            continue;

        if (!first)
            printf(" + ");

        fmpz_print(poly->coeffs + i);

        for (j = poly->nvars - 1; j >= 0; j--)
        {
            ulong m = poly->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            flint_printf("*%s^%wd", vars[dctx->perm[j]], e);
        }
        FLINT_ASSERT(k == 0);
        first = 0;
    }

    if (first)
        flint_printf("0");
}


void fmpz_mpolyd_print_simple(fmpz_mpolyd_t poly)
{
    int first = 0;
    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < poly->nvars; j++) {
        degb_prod *= poly->deg_bounds[j];
    }

    first = 1;
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (fmpz_is_zero(poly->coeffs + i))
            continue;

        if (!first)
            printf(" + ");

        fmpz_print(poly->coeffs + i);

        for (j = poly->nvars - 1; j >= 0; j--)
        {
            ulong m = poly->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            flint_printf("*x%d^%wd", j, e);
        }
        FLINT_ASSERT(k == 0);
        first = 0;
    }

    if (first)
        flint_printf("0");
}


void fmpz_mpoly_convert_from_fmpz_mpolyd(
                                  fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx,
                           const fmpz_mpolyd_t B, const fmpz_mpolyd_ctx_t dctx)
{
    slong i, j;
    slong degb_prod;
    slong * perm = dctx->perm;
    ulong * exps;
    TMP_INIT;

    FLINT_ASSERT(ctx->minfo->nvars == B->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < B->nvars; j++) {
        degb_prod *= B->deg_bounds[j];
    }

    TMP_START;
    exps = (ulong *) TMP_ALLOC(B->nvars*sizeof(ulong));

    fmpz_mpoly_zero(A, ctx);
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (fmpz_is_zero(B->coeffs + i))
            continue;

        for (j = B->nvars - 1; j >= 0; j--) 
        {
            ulong m = B->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            exps[perm[j]] = e;
        }
        FLINT_ASSERT(k == 0);

        fmpz_mpoly_set_coeff_fmpz_ui(A, B->coeffs + i, exps, ctx);
    }

    TMP_END;
}

int fmpz_mpolyd_CRT_nmod(fmpz_mpolyd_t A,
                                 const fmpz_mpolyd_t B, const fmpz_t Bm,
                                 const nmod_mpolyd_t C, const nmodf_ctx_t fctx)
{
    int Bok, Cok;
    slong carry;
    slong Bind, Cind;
    slong i, j;
    slong * inds;
    slong nvars = B->nvars;
    slong degb_prod;
    ulong hi, c;
    slong diff;
    fmpz_t zero;
    fmpz_t Bmn;
    slong * temp_deg_bounds;
    TMP_INIT;

    FLINT_ASSERT(B->nvars == C->nvars);

    TMP_START;
    temp_deg_bounds = (slong *) TMP_ALLOC(nvars*sizeof(slong));

    degb_prod = 1;
    diff = 0;
    for (j = 0; j < nvars; j++)
    {
        diff |= B->deg_bounds[j] - C->deg_bounds[j];
        temp_deg_bounds[j] = FLINT_MAX(B->deg_bounds[j], C->deg_bounds[j]);
        umul_ppmm(hi, degb_prod, degb_prod, temp_deg_bounds[j]);
        if (hi != WORD(0) || degb_prod < 0)
            return 0;
    }

    fmpz_init_set_ui(zero, 0);
    fmpz_init(Bmn);

    fmpz_mul_ui(Bmn, Bm, fctx->mod.n);
    c = fmpz_fdiv_ui(Bm, fctx->mod.n);
    c = n_invmod(c, fctx->mod.n);

    if (diff == 0) {
        /* both polynomials are packed into the same bounds */

        fmpz_mpolyd_set_nvars(A, nvars);
        fmpz_mpolyd_fit_length(A, degb_prod);
        for (j = 0; j < nvars; j++)
            A->deg_bounds[j] = temp_deg_bounds[j];

        for (i = 0; i < degb_prod; i++)
        {
            _fmpz_CRT_ui_precomp(A->coeffs + i, B->coeffs + i, Bm, C->coeffs[i], fctx->mod.n, fctx->mod.ninv, Bmn, c, 1);
        }

    } else {
        /* different bounds for packing */

        fmpz_mpolyd_t temp;
        fmpz_mpolyd_struct * T;

        if (A == B)
        {
            T = temp;
            fmpz_mpolyd_init(T, nvars);
        } else {
            T = A;
        }

        fmpz_mpolyd_set_nvars(T, nvars);
        fmpz_mpolyd_fit_length(T, degb_prod);
        for (j = 0; j < nvars; j++)
            T->deg_bounds[j] = temp_deg_bounds[j];

        inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
        for (j = 0; j < nvars; j++)
            inds[j] = 0;
        Bok = 1;
        Cok = 1;
        Bind = 0;
        Cind = 0;
        for (i = 0; i < degb_prod; i++)
        {
                   if (Bok && Cok) {
                _fmpz_CRT_ui_precomp(T->coeffs + i, B->coeffs + Bind++, Bm, C->coeffs[Cind++], fctx->mod.n, fctx->mod.ninv, Bmn, c, 1);
            } else if (Bok && !Cok) {
                _fmpz_CRT_ui_precomp(T->coeffs + i, B->coeffs + Bind++, Bm, 0                , fctx->mod.n, fctx->mod.ninv, Bmn, c, 1);
            } else if (!Bok && Cok) {
                _fmpz_CRT_ui_precomp(T->coeffs + i, zero              , Bm, C->coeffs[Cind++], fctx->mod.n, fctx->mod.ninv, Bmn, c, 1);
            } else {
                fmpz_zero(T->coeffs + i);
            }

            Bok = 1;
            Cok = 1;
            carry = 1;
            for (j = nvars - 1; j >= 0; j--)
            {
                inds[j] += carry;
                if (inds[j] < T->deg_bounds[j])
                {
                    carry = 0;
                    Bok = Bok && (inds[j] < B->deg_bounds[j]);
                    Cok = Cok && (inds[j] < C->deg_bounds[j]);
                } else
                {
                    carry = 1;
                    inds[j] = 0;
                }
            }
        }

        if (A == B)
        {
            fmpz_mpolyd_swap(A, T);
            fmpz_mpolyd_clear(T);
        } else {

        }
    }

    fmpz_clear(zero);
    fmpz_clear(Bmn);

    TMP_END;
    return 1;
}


void fmpz_mpolyd_height(fmpz_t max, fmpz_mpolyd_t A)
{
    slong degb_prod, i, j;
    fmpz_t t;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    fmpz_init(t);
    fmpz_zero(max);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_abs(t, A->coeffs + i);
        if (fmpz_cmp(max, t) < 0)
            fmpz_set(max, t);
    }

    fmpz_clear(t);
}

void fmpz_mpolyd_heights(fmpz_t max, fmpz_t sum, fmpz_mpolyd_t A)
{
    slong degb_prod, i, j;
    fmpz_t t;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    fmpz_init(t);
    fmpz_zero(max);
    fmpz_zero(sum);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_abs(t, A->coeffs + i);
        fmpz_add(sum, sum, t);
        if (fmpz_cmp(max, t) < 0)
            fmpz_set(max, t);
    }

    fmpz_clear(t);
}

void fmpz_mpolyd_divexact_fmpz_inplace(fmpz_mpolyd_t A, fmpz_t c)
{
    slong degb_prod, i, j;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    for (i = 0; i < degb_prod; i++)
    {
        fmpz_divexact(A->coeffs + i, A->coeffs + i, c);
    }
}

void fmpz_mpolyd_mul_scalar_inplace(fmpz_mpolyd_t A, fmpz_t c)
{
    slong degb_prod, i, j;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    for (i = 0; i < degb_prod; i++)
    {
        fmpz_mul(A->coeffs + i, A->coeffs + i, c);
    }
}

void fmpz_mpolyd_content(fmpz_t c, const fmpz_mpolyd_t A)
{
    slong degb_prod, j;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    _fmpz_vec_content(c, A->coeffs, degb_prod);
}

void fmpz_mpolyd_to_nmod_mpolyd(nmod_mpolyd_t Ap, fmpz_mpolyd_t A, const nmodf_ctx_t fctx)
{
    slong j;
    slong degb_prod;

    nmod_mpolyd_set_nvars(Ap, A->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        Ap->deg_bounds[j] = A->deg_bounds[j];
        degb_prod *= A->deg_bounds[j];
    }

    nmod_mpolyd_fit_length(Ap, degb_prod);
    _fmpz_vec_get_nmod_vec(Ap->coeffs, A->coeffs, degb_prod, fctx->mod);
}

void fmpz_mpolyd_set_nmod_mpolyd(fmpz_mpolyd_t A, nmod_mpolyd_t Ap, const nmodf_ctx_t fctx)
{
    slong j;
    slong degb_prod;

    fmpz_mpolyd_set_nvars(A, Ap->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < Ap->nvars; j++)
    {
        A->deg_bounds[j] = Ap->deg_bounds[j];
        degb_prod *= Ap->deg_bounds[j];
    }

    fmpz_mpolyd_fit_length(A, degb_prod);
    _fmpz_vec_set_nmod_vec(A->coeffs, Ap->coeffs, degb_prod, fctx->mod);
}



void fmpz_mpolyd_set(fmpz_mpolyd_t A, const fmpz_mpolyd_t B)
{
    slong j;
    slong degb_prod;

    fmpz_mpolyd_set_nvars(A, B->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < B->nvars; j++)
    {
        A->deg_bounds[j] = B->deg_bounds[j];
        degb_prod *= B->deg_bounds[j];
    }

    fmpz_mpolyd_fit_length(A, degb_prod);
    _fmpz_vec_set(A->coeffs, B->coeffs, degb_prod);
}


slong fmpz_mpolyd_leadmon(slong * exps, const fmpz_mpolyd_t A)
{

    slong i, j, k;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    for (i = degb_prod-1; i >= 0; i--)
    {
        if (!fmpz_is_zero(A->coeffs + i))
            break;
    }

    FLINT_ASSERT(i>=0);

    k = i;
    for (j = A->nvars - 1; j >= 0; j--) 
    {
        ulong m = A->deg_bounds[j];
        ulong e = k % m;
        k = k / m;
        exps[j] = e;
    }
    FLINT_ASSERT(k == 0);

    return i;
}

void fmpz_mpolyd_lc(fmpz_t a, const fmpz_mpolyd_t A)
{

    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    for (i = degb_prod-1; i >= 0; i--)
    {
        if (!fmpz_is_zero(A->coeffs + i))
        {
            fmpz_set(a, A->coeffs + i);
            return;
        }
    }
}



int fmpz_mpolyd_gcd_brown(fmpz_mpolyd_t G,
            fmpz_mpolyd_t Abar, fmpz_mpolyd_t Bbar,
                    fmpz_mpolyd_t A, fmpz_mpolyd_t B)
{
    int equal, success = 1;
    mp_limb_t p, old_p;
    slong j, nvars;
    slong lm_idx;
    slong * exp, * texp;
    fmpz_t gamma, m;
    fmpz_t gnm, gns, anm, ans, bnm, bns;
    fmpz_t lA, lB, cA, cB, cG, bound, temp, pp;
    nmod_mpolyd_t Gp, Apbar, Bpbar, Ap, Bp;
    nmodf_ctx_t fctx;
    TMP_INIT;

    TMP_START;

    nmodf_ctx_init(fctx, 2);
    nvars = A->nvars;

    nmod_mpolyd_init(Gp, nvars);
    nmod_mpolyd_init(Apbar, nvars);
    nmod_mpolyd_init(Bpbar, nvars);
    nmod_mpolyd_init(Ap, nvars);
    nmod_mpolyd_init(Bp, nvars);

    fmpz_init(cA);
    fmpz_init(cB);
    fmpz_init(cG);
    fmpz_init(lA);
    fmpz_init(lB);
    fmpz_init(gamma);
    fmpz_init(gnm);
    fmpz_init(gns);
    fmpz_init(anm);
    fmpz_init(ans);
    fmpz_init(bnm);
    fmpz_init(bns);
    fmpz_init(bound);
    fmpz_init(temp);
    fmpz_init_set_si(m, 1);
    fmpz_init(pp);

    fmpz_mpolyd_content(cA, A);
    fmpz_mpolyd_content(cB, B);
    fmpz_gcd(cG, cA, cB);
    fmpz_mpolyd_divexact_fmpz_inplace(A, cA);
    fmpz_mpolyd_divexact_fmpz_inplace(B, cB);
    fmpz_mpolyd_lc(lA, A);
    fmpz_mpolyd_lc(lB, B);
    fmpz_gcd(gamma, lA, lB);
    fmpz_mpolyd_height(bound, A);
    fmpz_mpolyd_height(temp, B);

    if (fmpz_cmp(bound, temp) < 0)
        fmpz_swap(bound, temp);
    fmpz_mul(bound, bound, gamma);
    fmpz_add(bound, bound, bound);

    exp = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    texp = (slong *) TMP_ALLOC(nvars*sizeof(slong));

    fmpz_mpolyd_leadmon(exp, A);
    fmpz_mpolyd_leadmon(texp, B);
    for (j = 0; j < nvars; j++)
        exp[j] = FLINT_MIN(exp[j], texp[j]);

    p = UWORD(1) << (FLINT_BITS - 1);

choose_next_prime:

    old_p = p;
    p = n_nextprime(p, 1);
    if (p <= old_p) {
        /* ran out of primes */
        success = 0;
        goto done;
    }
    fmpz_set_ui(pp, p);
    if (fmpz_divisible(lA, pp) || fmpz_divisible(lB, pp))
        goto choose_next_prime;

    nmodf_ctx_reset(fctx, p);
    fmpz_mpolyd_to_nmod_mpolyd(Ap, A, fctx);
    fmpz_mpolyd_to_nmod_mpolyd(Bp, B, fctx);
    success = nmod_mpolyd_gcd_brown_smprime(Gp, Apbar, Bpbar, Ap, Bp, fctx);
    if (!success)
        goto choose_next_prime;

    lm_idx = nmod_mpolyd_leadmon(texp, Gp);
    if (lm_idx <= 0)
    {
        /* Gp is 1, which means A and B are r.p. */
        FLINT_ASSERT(lm_idx == 0);
        fmpz_mpolyd_set_fmpz(G, cG);
        fmpz_mpolyd_swap(Abar, A);
        fmpz_divexact(temp, cA, cG);
        fmpz_mpolyd_mul_scalar_inplace(Abar, temp);
        fmpz_mpolyd_swap(Bbar, B);
        fmpz_divexact(temp, cB, cG);
        fmpz_mpolyd_mul_scalar_inplace(Bbar, temp);
        goto done;
    }

    equal = 1;
    for (j = 0; j < nvars; j++)
    {
        if (texp[j] > exp[j])
        {
            goto choose_next_prime;
        } else if (texp[j] < exp[j])
        {
            equal = 0;
            break;
        }
    }

    nmod_mpolyd_mul_scalar(Gp, fmpz_fdiv_ui(gamma, p), fctx);

    if (fmpz_is_one(m) || !equal)
    {
        fmpz_mpolyd_set_nmod_mpolyd(G, Gp, fctx);
        fmpz_mpolyd_set_nmod_mpolyd(Abar, Apbar, fctx);
        fmpz_mpolyd_set_nmod_mpolyd(Bbar, Bpbar, fctx);
        fmpz_set_ui(m, p);
        for (j = 0; j < nvars; j++)
            exp[j] = texp[j];

        goto choose_next_prime;
    }

    success = 1;
    success = success && fmpz_mpolyd_CRT_nmod(G, G, m, Gp, fctx);
    success = success && fmpz_mpolyd_CRT_nmod(Abar, Abar, m, Apbar, fctx);
    success = success && fmpz_mpolyd_CRT_nmod(Bbar, Bbar, m, Bpbar, fctx);
    fmpz_mul(m, m, pp);
    if (!success)
        goto done;

    if (fmpz_cmp(m, bound) <= 0)
        goto choose_next_prime;

    fmpz_mpolyd_heights(gnm, gns, G);
    fmpz_mpolyd_heights(anm, ans, Abar);
    fmpz_mpolyd_heights(bnm, bns, Bbar);
    fmpz_mul(ans, ans, gnm);
    fmpz_mul(anm, anm, gns);
    fmpz_mul(bns, bns, gnm);
    fmpz_mul(bnm, bnm, gns);

    if (fmpz_cmp(ans, anm) > 0)
        fmpz_swap(ans, anm);
    if (fmpz_cmp(bns, bnm) > 0)
        fmpz_swap(bns, bnm);
    fmpz_add(ans, ans, ans);
    fmpz_add(bns, bns, bns);
    if (fmpz_cmp(ans, m) >= 0 || fmpz_cmp(bns, m) >= 0)
        goto choose_next_prime;

    fmpz_mpolyd_content(temp, G);
    fmpz_mpolyd_divexact_fmpz_inplace(G, temp);
    fmpz_mpolyd_lc(temp, G);
    fmpz_mpolyd_divexact_fmpz_inplace(Abar, temp);
    fmpz_mpolyd_divexact_fmpz_inplace(Bbar, temp);

    fmpz_mpolyd_mul_scalar_inplace(G, cG);
    fmpz_divexact(temp, cA, cG);
    fmpz_mpolyd_mul_scalar_inplace(Abar, temp);
    fmpz_divexact(temp, cB, cG);
    fmpz_mpolyd_mul_scalar_inplace(Bbar, temp);

done:

    fmpz_clear(cA);
    fmpz_clear(cB);
    fmpz_clear(cG);
    fmpz_clear(lA);
    fmpz_clear(lB);
    fmpz_clear(gamma);
    fmpz_clear(gnm);
    fmpz_clear(gns);
    fmpz_clear(anm);
    fmpz_clear(ans);
    fmpz_clear(bnm);
    fmpz_clear(bns);
    fmpz_clear(bound);
    fmpz_clear(temp);
    fmpz_clear(m);
    fmpz_clear(pp);

    nmod_mpolyd_clear(Gp);
    nmod_mpolyd_clear(Apbar);
    nmod_mpolyd_clear(Bpbar);
    nmod_mpolyd_clear(Ap);
    nmod_mpolyd_clear(Bp);

    nmodf_ctx_clear(fctx);

    TMP_END;
    return success;
}


int fmpz_mpoly_gcd_brownold(fmpz_mpoly_t G,
                              const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mpolyd_t Ad, Bd, Gd, Abar, Bbar;
    fmpz_mpolyd_ctx_t dctx;
    slong nvars = ctx->minfo->nvars;

    success = 1;

    if (fmpz_mpoly_is_zero(A, ctx)) {
        fmpz_mpoly_set(G, B, ctx);
        goto cleanup_stage0;
    }
    if (fmpz_mpoly_is_zero(B, ctx)) {
        fmpz_mpoly_set(G, A, ctx);
        goto cleanup_stage0;
    }

    fmpz_mpolyd_ctx_init(dctx, nvars);
    success = fmpz_mpolyd_ctx_init_version1(dctx, A, B, ctx);
    if (!success)
    {
        fmpz_mpoly_zero(G, ctx);
        goto cleanup_stage1;
    }

    fmpz_mpolyd_init(Ad, nvars);
    fmpz_mpolyd_init(Bd, nvars);
    fmpz_mpolyd_init(Gd, nvars);
    fmpz_mpolyd_init(Abar, nvars);
    fmpz_mpolyd_init(Bbar, nvars);

    fmpz_mpoly_convert_to_fmpz_mpolyd(Ad, dctx, A, ctx);
    fmpz_mpoly_convert_to_fmpz_mpolyd(Bd, dctx, B, ctx);

    success = fmpz_mpolyd_gcd_brown(Gd, Abar, Bbar, Ad, Bd);
    if (!success)
    {
        fmpz_mpoly_zero(G, ctx);
    } else
    {
        fmpz_mpoly_convert_from_fmpz_mpolyd(G, ctx, Gd, dctx);
    }

    fmpz_mpolyd_clear(Bbar);
    fmpz_mpolyd_clear(Abar);
    fmpz_mpolyd_clear(Gd);
    fmpz_mpolyd_clear(Bd);
    fmpz_mpolyd_clear(Ad);

cleanup_stage1:

    fmpz_mpolyd_ctx_clear(dctx);

cleanup_stage0:

    if (success && (G->length > 0) && (fmpz_sgn(G->coeffs + 0) < 0))
        fmpz_mpoly_neg(G, G, ctx);

    return success;
}







void fmpz_mpolyu_height(
    fmpz_t max,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_t t;

    fmpz_init(t);
    fmpz_zero(max);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_abs(t, Ac->coeffs + j);
            if (fmpz_cmp(max, t) < 0)
                fmpz_set(max, t);
        }
    }

    fmpz_clear(t);
}


void fmpz_mpolyu_heights(
    fmpz_t max,
    fmpz_t sum,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_t t;

    fmpz_init(t);
    fmpz_zero(max);
    fmpz_zero(sum);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_abs(t, Ac->coeffs + j);
            fmpz_add(sum, sum, t);
            if (fmpz_cmp(max, t) < 0)
                fmpz_set(max, t);
        }
    }

    fmpz_clear(t);
}


/*
    A and B are assumed to be primitive and therefore are marked const
*/
int fmpz_mpolyu_gcd_brown(
    fmpz_mpolyu_t G,
    fmpz_mpolyu_t Abar,
    fmpz_mpolyu_t Bbar,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_t bound;
    slong offset, shift;
    mp_limb_t p, gammared;
    fmpz_t gamma, modulus;
    fmpz_t gnm, gns, anm, ans, bnm, bns;
    fmpz_t temp;
    fmpz_mpolyu_t T;
    nmod_mpolyun_t Gp, Abarp, Bbarp, Ap, Bp;
    nmod_mpoly_ctx_t pctx;
    mp_bitcnt_t bits = G->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, ctx->minfo->nvars - 1, G->bits, ctx->minfo);

    fmpz_init(gamma);
    fmpz_init(gnm);
    fmpz_init(gns);
    fmpz_init(anm);
    fmpz_init(ans);
    fmpz_init(bnm);
    fmpz_init(bns);
    fmpz_init(bound);
    fmpz_init(temp);
    fmpz_init_set_si(modulus, 1);

#if WANT_ASSERT
    fmpz_mpolyu_content_fmpz(temp, A, ctx);
    FLINT_ASSERT(fmpz_is_one(temp));
    fmpz_mpolyu_content_fmpz(temp, B, ctx);
    FLINT_ASSERT(fmpz_is_one(temp));
#endif

    fmpz_gcd(gamma, fmpz_mpolyu_leadcoeff_ref(A),
                    fmpz_mpolyu_leadcoeff_ref(B));

    fmpz_mpolyu_height(bound, A, ctx);
    fmpz_mpolyu_height(temp, B, ctx);
    if (fmpz_cmp(bound, temp) < 0)
        fmpz_swap(bound, temp);
    fmpz_mul(bound, bound, gamma);
    fmpz_add(bound, bound, bound);

    fmpz_mpolyu_init(T, bits, ctx);

    nmod_mpoly_ctx_init(pctx, ctx->minfo->nvars, ORD_LEX, 2);
    nmod_mpolyun_init(Ap, bits, pctx);
    nmod_mpolyun_init(Bp, bits, pctx);
    nmod_mpolyun_init(Gp, bits, pctx);
    nmod_mpolyun_init(Abarp, bits, pctx);
    nmod_mpolyun_init(Bbarp, bits, pctx);

    p = UWORD(1) << (FLINT_BITS - 1);

choose_prime:

    if (p >= UWORD_MAX_PRIME)
    {
        /* ran out of machine primes - absolute failure */
        success = 0;
        goto cleanup;
    }
    p = n_nextprime(p, 1);

    /* make sure reduction does not kill both lc(A) and lc(B) */
    gammared = fmpz_fdiv_ui(gamma, p);
    if (gammared == 0)
    {
        goto choose_prime;
    }

    nmod_mpoly_ctx_set_modulus(pctx, p);
    /* the unfortunate nmod poly's store their own context :( */
    nmod_mpolyun_set_mod(Ap, pctx->ffinfo->mod);
    nmod_mpolyun_set_mod(Bp, pctx->ffinfo->mod);
    nmod_mpolyun_set_mod(Gp, pctx->ffinfo->mod);
    nmod_mpolyun_set_mod(Abarp, pctx->ffinfo->mod);
    nmod_mpolyun_set_mod(Bbarp, pctx->ffinfo->mod);

    /* reduction should kill neither A nor B */
    fmpz_mpolyu_redto_nmod_mpolyun(Ap, pctx, A, ctx);
    fmpz_mpolyu_redto_nmod_mpolyun(Bp, pctx, B, ctx);
    FLINT_ASSERT(Ap->length > 0);
    FLINT_ASSERT(Bp->length > 0);

    success = nmod_mpolyun_gcd_brown_smprime(Gp, Abarp, Bbarp,
                                          Ap, Bp, ctx->minfo->nvars - 1, pctx);
    if (!success)
    {
        goto choose_prime;
    }

    if (nmod_mpolyun_is_nonzero_nmod(Gp, pctx))
    {
        fmpz_mpolyu_one(G, ctx);
        fmpz_mpolyu_set(Abar, A, ctx);
        fmpz_mpolyu_set(Bbar, B, ctx);
        goto successful_put_content;
    }

    if (!fmpz_is_one(modulus))
    {
        int cmp = 0;
        FLINT_ASSERT(G->length > 0);
        if (G->exps[0] != Gp->exps[0])
        {
            cmp = G->exps[0] > Gp->exps[0] ? 1 : -1;
        }
        if (cmp == 0)
        {
            slong k = nmod_poly_degree((Gp->coeffs + 0)->coeffs + 0);
            cmp = mpoly_monomial_cmp_nomask_extra(
                        (G->coeffs + 0)->exps + N*0,
                       (Gp->coeffs + 0)->exps + N*0, N, offset, k << shift);
        }

        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            fmpz_one(modulus);
        }
    }

    FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff(Gp, pctx));
    nmod_mpolyun_scalar_mul_nmod(Gp, gammared, pctx);

    if (!fmpz_is_one(modulus))
    {
        fmpz_mpolyu_addinterp_un(G, T, ctx, modulus, Gp, pctx);
        fmpz_mpolyu_addinterp_un(Abar, T, ctx, modulus, Abarp, pctx);
        fmpz_mpolyu_addinterp_un(Bbar, T, ctx, modulus, Bbarp, pctx);
    }
    else
    {
        fmpz_mpolyu_startinterp_un(G, ctx, Gp, pctx);
        fmpz_mpolyu_startinterp_un(Abar, ctx, Abarp, pctx);
        fmpz_mpolyu_startinterp_un(Bbar, ctx, Bbarp, pctx);
    }

    fmpz_mul_ui(modulus, modulus, p);

    if (fmpz_cmp(modulus, bound) <= 0)
    {
        goto choose_prime;
    }

    fmpz_mpolyu_heights(gnm, gns, G, ctx);
    fmpz_mpolyu_heights(anm, ans, Abar, ctx);
    fmpz_mpolyu_heights(bnm, bns, Bbar, ctx);
    fmpz_mul(ans, ans, gnm);
    fmpz_mul(anm, anm, gns);
    fmpz_mul(bns, bns, gnm);
    fmpz_mul(bnm, bnm, gns);

    if (fmpz_cmp(ans, anm) > 0)
        fmpz_swap(ans, anm);
    if (fmpz_cmp(bns, bnm) > 0)
        fmpz_swap(bns, bnm);
    fmpz_add(ans, ans, ans);
    fmpz_add(bns, bns, bns);
    if (fmpz_cmp(ans, modulus) < 0 && fmpz_cmp(bns, modulus) < 0)
    {
        goto successful;
    }

    /* do not reset modulus to 1 */
    goto choose_prime;

successful:

    FLINT_ASSERT(fmpz_equal(gamma, fmpz_mpolyu_leadcoeff_ref(G)));

    fmpz_mpolyu_content_fmpz(temp, G, ctx);
    fmpz_mpolyu_divexact_fmpz(G, G, temp, ctx);
    fmpz_mpolyu_divexact_fmpz(Abar, Abar, fmpz_mpolyu_leadcoeff_ref(G), ctx);
    fmpz_mpolyu_divexact_fmpz(Bbar, Bbar, fmpz_mpolyu_leadcoeff_ref(G), ctx);

successful_put_content:

    /* inputs were supposed to be primitive - nothing to do */
    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        fmpz_mul(temp, fmpz_mpolyu_leadcoeff_ref(G), fmpz_mpolyu_leadcoeff_ref(Abar));
        FLINT_ASSERT(fmpz_equal(temp, fmpz_mpolyu_leadcoeff_ref(A)));
        fmpz_mul(temp, fmpz_mpolyu_leadcoeff_ref(G), fmpz_mpolyu_leadcoeff_ref(Bbar));
        FLINT_ASSERT(fmpz_equal(temp, fmpz_mpolyu_leadcoeff_ref(B)));
    }
#endif

    fmpz_clear(gamma);
    fmpz_clear(gnm);
    fmpz_clear(gns);
    fmpz_clear(anm);
    fmpz_clear(ans);
    fmpz_clear(bnm);
    fmpz_clear(bns);
    fmpz_clear(bound);
    fmpz_clear(temp);
    fmpz_clear(modulus);

    nmod_mpolyun_clear(Gp, pctx);
    nmod_mpolyun_clear(Abarp, pctx);
    nmod_mpolyun_clear(Bbarp, pctx);
    nmod_mpolyun_clear(Ap, pctx);
    nmod_mpolyun_clear(Bp, pctx);

    nmod_mpoly_ctx_clear(pctx);

    fmpz_mpolyu_clear(T, ctx);

    return success;
}





int fmpz_mpoly_gcd_brown(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong * perm;
    ulong * shift, * stride;
    slong i;
    mp_bitcnt_t new_bits;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    fmpz_t cA, cB, cG;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }
        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars == 1)
    {
        slong shiftA, shiftB;
        fmpz_poly_t a, b, g;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        fmpz_mpoly_to_fmpz_poly_keepbits(a, &shiftA, A, 0, ctx);
        fmpz_mpoly_to_fmpz_poly_keepbits(b, &shiftB, B, 0, ctx);
        fmpz_poly_gcd(g, a, b);
        fmpz_mpoly_from_fmpz_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);
        return 1;
    }

    perm = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i + 1 < ctx->minfo->nvars ? i + 1 : 0;
        shift[i] = 0;
        stride[i] = 1;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX);
    fmpz_mpolyu_init(Au, new_bits, uctx);
    fmpz_mpolyu_init(Bu, new_bits, uctx);
    fmpz_mpolyu_init(Gu, new_bits, uctx);
    fmpz_mpolyu_init(Abaru, new_bits, uctx);
    fmpz_mpolyu_init(Bbaru, new_bits, uctx);

    fmpz_init(cA);
    fmpz_init(cB);
    fmpz_init(cG);

    fmpz_mpoly_to_mpolyu_perm_deflate(Au, A, perm, shift, stride, uctx, ctx);
    fmpz_mpoly_to_mpolyu_perm_deflate(Bu, B, perm, shift, stride, uctx, ctx);

    fmpz_mpolyu_content_fmpz(cA, Au, uctx);
    fmpz_mpolyu_content_fmpz(cB, Bu, uctx);
    fmpz_gcd(cG, cA, cB);
    fmpz_mpolyu_divexact_fmpz(Au, Au, cA, uctx);
    fmpz_mpolyu_divexact_fmpz(Bu, Bu, cB, uctx);
    success = fmpz_mpolyu_gcd_brown(Gu, Abaru, Bbaru, Au, Bu, uctx);

    if (success)
    {
        fmpz_mpoly_from_mpolyu_perm_inflate(G, new_bits, Gu, perm, shift, stride, uctx, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
            fmpz_neg(cG, cG);
        fmpz_mpoly_scalar_mul_fmpz(G, G, cG, ctx);
    }

    fmpz_clear(cA);
    fmpz_clear(cB);
    fmpz_clear(cG);

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Gu, uctx);
    fmpz_mpolyu_clear(Abaru, uctx);
    fmpz_mpolyu_clear(Bbaru, uctx);
    fmpz_mpoly_ctx_clear(uctx);

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    return success;
}
