/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


void nmod_mpolyd_ctx_init(nmod_mpolyd_ctx_t dctx, slong nvars)
{
    slong i;

    dctx->nvars = nvars;
    dctx->perm = (slong *) flint_malloc(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }
}


int nmod_mpolyd_ctx_settle(nmod_mpolyd_ctx_t dctx,
                            const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, degb_prod;
    slong * Aexps, * Bexps, * deg_bounds;
    slong nvars = ctx->minfo->nvars;
    slong * perm = dctx->perm;
    TMP_INIT;

    success = 0;

    TMP_START;
    Aexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Bexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
        goto cleanup;
    nmod_mpoly_degrees_si(Aexps, A, ctx);
    nmod_mpoly_degrees_si(Bexps, B, ctx);

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

void nmod_mpolyd_ctx_clear(nmod_mpolyd_ctx_t dctx)
{
    flint_free(dctx->perm);
}

void nmod_mpolyd_init(nmod_mpolyd_t poly, slong nvars)
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
    poly->coeffs = (mp_limb_t *) flint_malloc(poly->coeff_alloc*sizeof(mp_limb_t));
    for (i = 0; i < poly->coeff_alloc; i++)
    {
        poly->coeffs[i] = UWORD(0);
    }
}

void nmod_mpolyd_fit_length(nmod_mpolyd_t poly, slong len) {
    if (poly->coeff_alloc < len) {
        poly->coeffs = (mp_limb_t *) flint_realloc(poly->coeffs, len*sizeof(mp_limb_t));
        poly->coeff_alloc = len;
    }
}

void nmod_mpolyd_set_nvars(nmod_mpolyd_t poly, slong nvars) {

    poly->nvars = nvars;
    if (poly->degb_alloc < nvars) {
        poly->deg_bounds = (slong *) flint_realloc(poly->deg_bounds, nvars*sizeof(slong));
        poly->degb_alloc = nvars;
    }
}

void nmod_mpolyd_zero(nmod_mpolyd_t poly)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeffs[0] = UWORD(0);
}



void nmod_mpolyd_clear(nmod_mpolyd_t poly)
{
    flint_free(poly->deg_bounds);
    flint_free(poly->coeffs);
    poly->deg_bounds = NULL;
    poly->coeffs = NULL;
}

void nmod_mpoly_convert_to_nmod_mpolyd(
                            nmod_mpolyd_t poly1, const nmod_mpolyd_ctx_t dctx,
                          const nmod_mpoly_t poly2, const nmod_mpoly_ctx_t ctx)
{
    slong degb_prod;
    slong i, j, N;
    slong * exps;
    const slong * perm = dctx->perm;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    nmod_mpolyd_set_nvars(poly1, ctx->minfo->nvars);

    FLINT_ASSERT(poly2->bits <= FLINT_BITS)

    if (poly2->length == 0)
    {
        nmod_mpolyd_zero(poly1);
        return;
    }

    TMP_START;
    exps = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));

    nmod_mpoly_degrees_si(exps, poly2, ctx);
    degb_prod = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        poly1->deg_bounds[i] = exps[perm[i]] + 1;
        degb_prod *= poly1->deg_bounds[i];
    }

    nmod_mpolyd_fit_length(poly1, degb_prod);
    for (i = 0; i < degb_prod; i++)
    {
        poly1->coeffs[i] = UWORD(0);
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

        poly1->coeffs[off] = poly2->coeffs[i];
    }

    TMP_END;
}

void nmod_mpolyd_print(nmod_mpolyd_t poly, const char ** vars,
                                                  const nmod_mpolyd_ctx_t dctx)
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

        if (poly->coeffs[i] == 0)
            continue;

        if (!first)
            printf(" + ");

        flint_printf("%wu", poly->coeffs[i]);

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

void nmod_mpolyd_print_simple(nmod_mpolyd_t poly)
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

        if (poly->coeffs[i] == 0)
            continue;

        if (!first)
            printf(" + ");

        flint_printf("%wu", poly->coeffs[i]);

        for (j = poly->nvars - 1; j >= 0; j--) 
        {
            ulong m = poly->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            flint_printf("*x%wd^%wd", j, e);
        }
        FLINT_ASSERT(k == 0);
        first = 0;
    }

    if (first)
        flint_printf("0");
}

void nmod_mpolyd_add(nmod_mpolyd_t A, const nmod_mpolyd_t B,
                                 const nmod_mpolyd_t C, const nmodf_ctx_t fctx)
{
    int Bok, Cok;
    slong Bind, Cind;
    slong i, j;
    slong * inds;
    slong nvars = B->nvars;
    slong degb_prod;
    TMP_INIT;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);
    FLINT_ASSERT(B->nvars == C->nvars);

    nmod_mpolyd_set_nvars(A, B->nvars);

    degb_prod = 1;
    for (j = 0; j < nvars; j++)
    {
        A->deg_bounds[j] = FLINT_MAX(B->deg_bounds[j], C->deg_bounds[j]);
        degb_prod *= A->deg_bounds[j];
    }
    nmod_mpolyd_fit_length(A, degb_prod);

    TMP_START;
    inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (j = 0; j < nvars; j++)
        inds[j] = 0;
    Bok = 1;
    Cok = 1;
    Bind = 0;
    Cind = 0;
    for (i = 0; i < degb_prod; i++)
    {
        mp_limb_t x = WORD(0);
        int carry = 1;

        if (Bok)
            x += B->coeffs[Bind++];
        if (Cok)
            x = nmod_add(x, C->coeffs[Cind++], fctx->mod);

        A->coeffs[i] = x;

        Bok = 1;
        Cok = 1;
        for (j = nvars - 1; j >= 0; j--)
        {
            inds[j] += carry;
            if (inds[j] < A->deg_bounds[j])
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

    TMP_END;
}


void nmod_mpolyd_sub(nmod_mpolyd_t A, const nmod_mpolyd_t B,
                                 const nmod_mpolyd_t C, const nmodf_ctx_t fctx)
{
    int Bok, Cok;
    slong Bind, Cind;
    slong i, j;
    slong * inds;
    slong nvars = B->nvars;
    slong degb_prod;
    TMP_INIT;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);
    FLINT_ASSERT(B->nvars == C->nvars);

    nmod_mpolyd_set_nvars(A, B->nvars);

    degb_prod = 1;
    for (j = 0; j < nvars; j++)
    {
        A->deg_bounds[j] = FLINT_MAX(B->deg_bounds[j], C->deg_bounds[j]);
        degb_prod *= A->deg_bounds[j];
    }
    nmod_mpolyd_fit_length(A, degb_prod);

    TMP_START;
    inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (j = 0; j < nvars; j++)
        inds[j] = 0;
    Bok = 1;
    Cok = 1;
    Bind = 0;
    Cind = 0;
    for (i = 0; i < degb_prod; i++)
    {
        mp_limb_t x = WORD(0);
        int carry = 1;

        if (Bok)
            x += B->coeffs[Bind++];
        if (Cok)
            x = nmod_sub(x, C->coeffs[Cind++], fctx->mod);

        A->coeffs[i] = x;

        Bok = 1;
        Cok = 1;
        for (j = nvars - 1; j >= 0; j--)
        {
            inds[j] += carry;
            if (inds[j] < A->deg_bounds[j])
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

    TMP_END;
}




void nmod_mpoly_convert_from_nmod_mpolyd(
                                  nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx,
                           const nmod_mpolyd_t B, const nmod_mpolyd_ctx_t dctx)
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

    nmod_mpoly_zero(A, ctx);
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (B->coeffs[i] == 0)
            continue;

        for (j = B->nvars - 1; j >= 0; j--) 
        {
            ulong m = B->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            exps[perm[j]] = e;
        }
        FLINT_ASSERT(k == 0);

        nmod_mpoly_set_term_ui_ui(A, B->coeffs[i], exps, ctx);
    }

    TMP_END;
}






void nmod_mpolyd_gcd_brown_univar(nmod_mpolyd_t G,
     nmod_mpolyd_t Abar, nmod_mpolyd_t Bbar,
          const nmod_mpolyd_t A, const nmod_mpolyd_t B,
                 const nmodf_ctx_t fctx)
{
        slong Alen, Blen;

        Alen = A->deg_bounds[0];
        while (A->coeffs[Alen-1] == 0) {
            Alen --;
            if (Alen == 0)
                break;
        }

        Blen = B->deg_bounds[0];
        while (B->coeffs[Blen-1] == 0) {
            Blen --;
            if (Blen == 0)
                break;
        }

        nmod_mpolyd_set_nvars(G, 1);
        nmod_mpolyd_set_nvars(Abar, 1);
        nmod_mpolyd_set_nvars(Bbar, 1);
        if (Alen == 0)
        {
            if (Blen == 0)
            {
                G->coeffs[0] = 0;
                G->deg_bounds[0] = 1;
                Abar->coeffs[0] = 0;
                Abar->deg_bounds[0] = 1;
                Bbar->coeffs[0] = 0;
                Bbar->deg_bounds[0] = 1;
            } else
            {
                nmod_mpolyd_fit_length(G, Blen);
                _nmod_poly_make_monic(G->coeffs, B->coeffs, Blen, fctx->mod);
                G->deg_bounds[0] = Blen;
                Abar->coeffs[0] = 0;
                Abar->deg_bounds[0] = 1;
                Bbar->coeffs[0] = 1;
                Bbar->deg_bounds[0] = 1;
            }

        } else {
            if (Blen == 0) {
                nmod_mpolyd_fit_length(G, Alen);
                _nmod_poly_make_monic(G->coeffs, A->coeffs, Alen, fctx->mod);
                G->deg_bounds[0] = Alen;
                Abar->coeffs[0] = 1;
                Abar->deg_bounds[0] = 1;
                Bbar->coeffs[0] = 0;
                Bbar->deg_bounds[0] = 1;
            } else {
                if (Alen >= Blen) {
                    nmod_mpolyd_fit_length(G, Blen);
                    G->deg_bounds[0] = _nmod_poly_gcd_euclidean(G->coeffs,
                                             A->coeffs, Alen, B->coeffs, Blen,
                                                                    fctx->mod);
                } else {
                    nmod_mpolyd_fit_length(G, Alen);
                    G->deg_bounds[0] = _nmod_poly_gcd_euclidean(G->coeffs,
                                             B->coeffs, Blen, A->coeffs, Alen,
                                                                    fctx->mod);
                }
                if (G->deg_bounds[0] <= 1)
                {
                    G->coeffs[0] = 1;
                } else
                {
                    _nmod_poly_make_monic(G->coeffs,
                                       G->coeffs, G->deg_bounds[0], fctx->mod);
                }

                Abar->deg_bounds[0] = Alen - G->deg_bounds[0] + 1;
                nmod_mpolyd_fit_length(Abar, Alen - G->deg_bounds[0] + 1);
                _nmod_poly_div(Abar->coeffs, A->coeffs, Alen,
                                       G->coeffs, G->deg_bounds[0], fctx->mod);

                Bbar->deg_bounds[0] = Blen - G->deg_bounds[0] + 1;
                nmod_mpolyd_fit_length(Bbar, Blen - G->deg_bounds[0] + 1);
                _nmod_poly_div(Bbar->coeffs, B->coeffs, Blen,
                                       G->coeffs, G->deg_bounds[0], fctx->mod);
            }
        }
        return;
}



void nmod_mpolyd_last_content(nmod_poly_t cont, const nmod_mpolyd_t A,
                                                        const nmodf_ctx_t fctx) 
{
    mp_limb_t * temp;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
    slong i, j, Plen;
    slong degb_prod, degb_last=0;
    TMP_INIT;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }

    TMP_START;
    nmod_poly_zero(cont);
    nmod_poly_fit_length(cont, degb_last);
    temp = (mp_limb_t *) TMP_ALLOC(degb_last*sizeof(mp_limb_t));

    for (i = 0; i < degb_prod; i += degb_last)
    {
        mp_limb_t * P = A->coeffs + i;
        Plen = degb_last;
        while (P[Plen-1] == 0)
        {
            Plen --;
            if (Plen == 0)
                break;
        }

        if (Plen != 0)
        {
            if (cont->length == 0)
            {
                _nmod_poly_make_monic(cont->coeffs, P, Plen, fctx->mod);
                cont->length = Plen;
            } else
            {
                if (cont->length < Plen) {
                    cont->length = _nmod_poly_gcd_euclidean(temp,
                                         P, Plen, cont->coeffs, cont->length,
                                                                    fctx->mod);
                } else {
                    cont->length = _nmod_poly_gcd_euclidean(temp,
                                         cont->coeffs, cont->length, P, Plen,
                                                                    fctx->mod);
                }
                if (cont->length <= 1)
                {
                    cont->coeffs[0] = 1;
                } else
                {
                    _nmod_poly_make_monic(cont->coeffs, temp, cont->length, fctx->mod);
                }
            }
        }
    }
    TMP_END;
}

slong nmod_mpolyd_last_degree(const nmod_mpolyd_t A, const nmodf_ctx_t fctx) 
{                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    slong i, j, Plen, degree;
    slong degb_prod, degb_last=0;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }

    degree = -1;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        mp_limb_t * P = A->coeffs + i;
        Plen = degb_last;
        while (P[Plen-1] == 0)
        {
            Plen --;
            if (Plen == 0)
                break;
        }
        degree = FLINT_MAX(degree, Plen - 1);
    }
    return degree;
}


void nmod_mpolyd_div_last_poly(nmod_mpolyd_t A, nmod_poly_t b,
                                                        const nmodf_ctx_t fctx)
{
    slong i, j, k;
    slong degb_prod, degb_last=0, new_degb_last;
    mp_limb_t * temp;
    TMP_INIT;

    FLINT_ASSERT(b->length != 0);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }
    new_degb_last = degb_last - (b->length - 1);

    TMP_START;
    temp = (mp_limb_t *) TMP_ALLOC(new_degb_last*sizeof(mp_limb_t));

    k = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        mp_limb_t * P = A->coeffs + i;
        slong Plen = degb_last;
        while (P[Plen-1] == 0)
        {
            Plen --;
            if (Plen == 0)
                break;
        }

        j = 0;
        if (Plen != 0)
        {
            FLINT_ASSERT(Plen > b->length - 1);
            _nmod_poly_div(temp, P, Plen, b->coeffs, b->length, fctx->mod);
            while (j < Plen - (b->length - 1))
            {
                A->coeffs[k++] = temp[j];
                j++;
            }
        }
        while (j < new_degb_last)
        {
            A->coeffs[k++] = WORD(0);
            j++;
        }
    }

    A->deg_bounds[A->nvars - 1] = new_degb_last;

    TMP_END;
}

void nmod_mpolyd_mul_last_poly(nmod_mpolyd_t M,
                   nmod_mpolyd_t A, nmod_poly_t b, const nmodf_ctx_t fctx)
{
    slong i, j, k;
    slong degb_prod, degb_prod_last=0, degb_last=0, new_degb_last;

    FLINT_ASSERT(M != A);
    FLINT_ASSERT(b->length != 0);

    nmod_mpolyd_set_nvars(M, A->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        M->deg_bounds[j] = A->deg_bounds[j];
        degb_last = A->deg_bounds[j];
        degb_prod_last = degb_prod;
        degb_prod *= degb_last;
    }
    new_degb_last = degb_last + (b->length - 1);

    M->deg_bounds[M->nvars - 1] = degb_last + (b->length - 1);
    nmod_mpolyd_fit_length(M, degb_prod_last*(degb_last + (b->length - 1)));

    k = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        mp_limb_t * P = A->coeffs + i;
        slong Plen = degb_last;
        while (P[Plen-1] == 0)
        {
            Plen --;
            if (Plen == 0)
                break;
        }

        j = 0;
        if (Plen != 0)
        {
            _nmod_poly_mul_classical(M->coeffs + k, P, Plen,
                                              b->coeffs, b->length, fctx->mod);
            j +=  Plen + (b->length - 1);
        }
        while (j < new_degb_last)
        {

            M->coeffs[k+j] = WORD(0);
            j++;
        }

        k += new_degb_last;
    }
}


void nmod_mpolyd_mul_scalar(nmod_mpolyd_t A, mp_limb_t b, const nmodf_ctx_t fctx)
{
    slong j, degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_prod *= A->deg_bounds[j];
    }

    _nmod_vec_scalar_mul_nmod(A->coeffs, A->coeffs, degb_prod, b, fctx->mod);
}


void nmod_mpolyd_last_lc(nmod_poly_t lc, const nmod_mpolyd_t A,
                                                        const nmodf_ctx_t fctx)
{

    slong i, j, lc_len;
    slong degb_prod, degb_last=0;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }

    nmod_poly_zero(lc);

    for (i = degb_prod - degb_last; i >= 0; i -= degb_last)
    {
        lc_len = degb_last;
        while (lc_len > 0)
        {
            if (A->coeffs[i+lc_len-1] != 0)
            {
                nmod_poly_fit_length(lc, lc_len);
                flint_mpn_copyi(lc->coeffs, A->coeffs + i, lc_len);
                lc->length = lc_len;
                return;
            }
            lc_len--;
        }
    }

    return;
}


void nmod_mpolyd_eval_last(nmod_mpolyd_t E, const nmod_mpolyd_t A,
                                      mp_limb_t alpha, const nmodf_ctx_t fctx)
{

    slong i, j, k;
    slong degb_prod, degb_last=0, degb_prod_last=0;

    nmod_mpolyd_set_nvars(E, A->nvars - 1);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_last = A->deg_bounds[j];
        degb_prod_last = degb_prod;
        degb_prod *= degb_last;
        if (j < E->nvars)
        {
            E->deg_bounds[j] = A->deg_bounds[j];
        }
    }

    nmod_mpolyd_fit_length(E, degb_prod_last);

    j = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        mp_limb_t pp0, pp1, v = WORD(0);
        k = degb_last;
        while (--k >= 0)
        {
            umul_ppmm(pp1, pp0, v, alpha);
            add_ssaaaa(pp1, pp0, pp1, pp0, WORD(0), A->coeffs[i + k]);
            NMOD_RED2(v, pp1, pp0, fctx->mod);
        }
        E->coeffs[j++] = v;
    }

    FLINT_ASSERT(j == degb_prod_last);

    return;
}

slong nmod_mpolyd_leadmon(slong * exps, const nmod_mpolyd_t A)
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
        if (A->coeffs[i] != 0)
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


void nmod_mpolyd_set_ui(nmod_mpolyd_t A, mp_limb_t v)
{
    slong j;

    for (j = 0; j < A->nvars; j++)
        A->deg_bounds[j] = WORD(1);

    A->coeffs[0] = v;
}

void nmod_mpolyd_set(nmod_mpolyd_t A, const nmod_mpolyd_t B)
{
    slong j;
    slong degb_prod;

    nmod_mpolyd_set_nvars(A, B->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < B->nvars; j++)
    {
        A->deg_bounds[j] = B->deg_bounds[j];
        degb_prod *= B->deg_bounds[j];
    }

    nmod_mpolyd_fit_length(A, degb_prod);
    flint_mpn_copyi(A->coeffs, B->coeffs, degb_prod);
}


void nmod_mpolyd_interpolation(nmod_mpolyd_t fxn,
             nmod_mpolyd_t newvalue, nmod_poly_t modulus, mp_limb_t alpha,
                 const nmodf_ctx_t fctx)
{
    nmod_mpolyd_t temp;
    nmod_mpolyd_t E;
    slong nvars = fxn->nvars;

    nmod_mpolyd_init(E, nvars);
    nmod_mpolyd_eval_last(E, fxn, alpha, fctx);

    nmod_mpolyd_init(temp, nvars);
    nmod_mpolyd_sub(temp, newvalue, E, fctx);

    FLINT_ASSERT(temp->nvars = nvars - 1);
    nmod_mpolyd_set_nvars(temp, nvars);
    temp->deg_bounds[nvars - 1] = 1;

    nmod_mpolyd_mul_last_poly(E, temp, modulus, fctx);
    nmod_mpolyd_add(temp, fxn, E, fctx);
    nmod_mpolyd_set(fxn, temp);
    nmod_mpolyd_clear(temp);
    nmod_mpolyd_clear(E);
}

void nmod_mpolyd_promote(nmod_mpolyd_t fxn,
             nmod_mpolyd_t newvalue)
{
    slong nvars = newvalue->nvars;

    nmod_mpolyd_set_nvars(newvalue, nvars+1);
    newvalue->deg_bounds[nvars] = 1;

    nmod_mpolyd_set(fxn, newvalue);
}



int nmod_mpolyd_gcd_brown(nmod_mpolyd_t G,
     nmod_mpolyd_t Abar, nmod_mpolyd_t Bbar,
            nmod_mpolyd_t A, nmod_mpolyd_t B,
                 const nmodf_ctx_t fctx) {

    int success;
    slong i, j, bound;
    slong nvars = A->nvars;
    mp_limb_t alpha, gamma_eval;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, lcA, lcB, gamma;
    nmod_poly_t cGs, cAbars, cBbars, modulus, modulus2;
    nmod_mpolyd_t Gs, Abars, Bbars, phiA, phiB, gs, abars, bbars;
    slong leadmon_gs_idx;
    slong * leadmon_gs, * leadmon_Gs;
    slong deggamma, degGs, degA, degB, degAbars, degBbars;

    FLINT_ASSERT(G != A);
    FLINT_ASSERT(G != B);
    FLINT_ASSERT(A->nvars == B->nvars);

    if (A->nvars == 1) {
        nmod_mpolyd_gcd_brown_univar(G, Abar, Bbar, A, B, fctx);
        return 1;
    }

    leadmon_gs = (slong *) flint_malloc(nvars*sizeof(slong));
    leadmon_Gs = (slong *) flint_malloc(nvars*sizeof(slong));

    nmod_poly_init(cA, fctx->mod.n);
    nmod_poly_init(cB, fctx->mod.n);
    nmod_mpolyd_last_content(cA, A, fctx);
    nmod_mpolyd_last_content(cB, B, fctx);

    nmod_mpolyd_div_last_poly(A, cA, fctx);
    nmod_mpolyd_div_last_poly(B, cB, fctx);

    nmod_poly_init(cG, fctx->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, fctx->mod.n);
    nmod_poly_init(cBbar, fctx->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(lcA, fctx->mod.n);
    nmod_poly_init(lcB, fctx->mod.n);
    nmod_mpolyd_last_lc(lcA, A, fctx);
    nmod_mpolyd_last_lc(lcB, B, fctx);

    nmod_poly_init(gamma, fctx->mod.n);
    nmod_poly_gcd_euclidean(gamma, lcA, lcB);

    bound = 1 + nmod_poly_degree(gamma)
               + FLINT_MAX(nmod_mpolyd_last_degree(A, fctx),
                           nmod_mpolyd_last_degree(B, fctx));

    nmod_mpolyd_init(Gs, nvars);
    nmod_mpolyd_init(Abars, nvars);
    nmod_mpolyd_init(Bbars, nvars);
    nmod_mpolyd_init(phiA, nvars - 1);
    nmod_mpolyd_init(phiB, nvars - 1);

    nmod_mpolyd_init(gs, nvars - 1);
    nmod_mpolyd_init(abars, nvars - 1);
    nmod_mpolyd_init(bbars, nvars - 1);

    nmod_poly_init(modulus, fctx->mod.n);
    nmod_poly_init(modulus2, fctx->mod.n);

    i = 0;

    for (alpha = 0; alpha < fctx->mod.n; alpha++)
    {
        gamma_eval = nmod_poly_evaluate_nmod(gamma, alpha);
        if (gamma_eval == 0)
            goto break_continue;

        nmod_mpolyd_eval_last(phiA, A, alpha, fctx);
        nmod_mpolyd_eval_last(phiB, B, alpha, fctx);

        success = nmod_mpolyd_gcd_brown(gs, abars, bbars, phiA, phiB, fctx);
        if (success == 0)
            goto break_continue;

        leadmon_gs_idx = nmod_mpolyd_leadmon(leadmon_gs, gs);

        if (leadmon_gs_idx <= 0)
        {
            FLINT_ASSERT(leadmon_gs_idx == 0);
            nmod_mpolyd_set_ui(Gs, WORD(1));
            nmod_mpolyd_set(Abars, A);
            nmod_mpolyd_set(Bbars, B);
            goto successful;        
        }

        if (i > 0)
        {
            nmod_mpolyd_leadmon(leadmon_Gs, Gs);
            for (j = 0; j < nvars - 1; j++)
            {
                if (leadmon_gs[j] > leadmon_Gs[j])
                {
                    goto break_continue;
                } else if (leadmon_gs[j] < leadmon_Gs[j])
                {
                    nmod_mpolyd_zero(Gs);
                    nmod_mpolyd_zero(Abars);
                    nmod_mpolyd_zero(Bbars);
                    i = 0;
                }
            }
        }

        nmod_mpolyd_mul_scalar(gs, gamma_eval, fctx);

        if (i > 0)
        {
            nmod_poly_scalar_mul_nmod(modulus, modulus,
                n_invmod(nmod_poly_evaluate_nmod(modulus, alpha), fctx->mod.n));
            nmod_mpolyd_interpolation(Gs, gs, modulus, alpha, fctx);
            nmod_mpolyd_interpolation(Abars, abars, modulus, alpha, fctx);
            nmod_mpolyd_interpolation(Bbars, bbars, modulus, alpha, fctx);

        } else
        {
            nmod_poly_one(modulus);
            nmod_mpolyd_promote(Gs, gs);
            nmod_mpolyd_promote(Abars, abars);
            nmod_mpolyd_promote(Bbars, bbars);
        }

        nmod_poly_scalar_mul_nmod(modulus2, modulus, alpha);
        nmod_poly_shift_left(modulus, modulus, 1);
        nmod_poly_sub(modulus, modulus, modulus2);

        i++;
        if (i < bound)
            continue;

        deggamma = nmod_poly_degree(gamma);
        degGs = nmod_mpolyd_last_degree(Gs, fctx);
        degA = nmod_mpolyd_last_degree(A, fctx);
        degB = nmod_mpolyd_last_degree(B, fctx);
        degAbars = nmod_mpolyd_last_degree(Abars, fctx);
        degBbars = nmod_mpolyd_last_degree(Bbars, fctx);

        if (   deggamma + degA == degGs + degAbars
            && deggamma + degB == degGs + degBbars
           )
        {
            goto successful;
        } else
        {
            nmod_mpolyd_zero(Gs);
            nmod_mpolyd_zero(Abars);
            nmod_mpolyd_zero(Bbars);
            i = 0;
            continue;
        }

break_continue:
        (void)(NULL);
    }

    success = 0;

cleanup:

    flint_free(leadmon_gs);
    flint_free(leadmon_Gs);

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);

    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(lcA);
    nmod_poly_clear(lcB);

    nmod_poly_clear(gamma);

    nmod_mpolyd_clear(Gs);
    nmod_mpolyd_clear(Abars);
    nmod_mpolyd_clear(Bbars);
    nmod_mpolyd_clear(phiA);
    nmod_mpolyd_clear(phiB);

    nmod_mpolyd_clear(gs);
    nmod_mpolyd_clear(abars);
    nmod_mpolyd_clear(bbars);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;

successful:

    nmod_poly_init(cGs, fctx->mod.n);
    nmod_poly_init(cAbars, fctx->mod.n);
    nmod_poly_init(cBbars, fctx->mod.n);
    nmod_mpolyd_last_content(cGs, Gs, fctx);
    nmod_mpolyd_last_content(cAbars, Abars, fctx);
    nmod_mpolyd_last_content(cBbars, Bbars, fctx);

    nmod_mpolyd_div_last_poly(Gs, cGs, fctx);
    nmod_mpolyd_div_last_poly(Abars, cAbars, fctx);
    nmod_mpolyd_div_last_poly(Bbars, cBbars, fctx);

    nmod_mpolyd_mul_last_poly(G, Gs, cG, fctx);
    nmod_mpolyd_mul_last_poly(Abar, Abars, cAbar, fctx);
    nmod_mpolyd_mul_last_poly(Bbar, Bbars, cBbar, fctx);

    nmod_poly_clear(cGs);
    nmod_poly_clear(cAbars);
    nmod_poly_clear(cBbars);

    success = 1;
    goto cleanup;
}


