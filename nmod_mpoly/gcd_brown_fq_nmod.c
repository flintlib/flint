/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_poly.h"

/*
    helper functions for
    small prime version of the modular gcd over finite fields
*/

void fq_nmod_mpolyd_ctx_init(fq_nmod_mpolyd_ctx_t dctx, slong nvars,
                                                        mp_limb_t p, slong deg)
{
    slong i;
    fmpz_t P;
    fmpz_init_set_ui(P, p);

    dctx->nvars = nvars;
    dctx->perm = (slong *) flint_malloc(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }

    fq_nmod_ctx_init(dctx->fqctx, P, deg, "#");
    fmpz_clear(P);
}


void fq_nmod_mpolyd_ctx_clear(fq_nmod_mpolyd_ctx_t dctx)
{
    flint_free(dctx->perm);
    fq_nmod_ctx_clear(dctx->fqctx);
}

void fq_nmod_mpolyd_init(fq_nmod_mpolyd_t poly, slong nvars,
                                               const fq_nmod_mpolyd_ctx_t dctx)
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
    poly->coeffs = (fq_nmod_struct *) flint_malloc(poly->coeff_alloc
                                                      *sizeof(fq_nmod_struct));
    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fq_nmod_init(poly->coeffs + i, dctx->fqctx);
    }
}

NMOD_MPOLY_INLINE void fq_nmod_mpolyd_swap(fq_nmod_mpolyd_t poly1,
                                fq_nmod_mpolyd_t poly2)
{
   fq_nmod_mpolyd_struct t = *poly1;
   *poly1 = *poly2;
   *poly2 = t;
}

void fq_nmod_mpolyd_fit_length(fq_nmod_mpolyd_t poly, slong len,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{
    slong i;
    if (poly->coeff_alloc < len) {
        slong new_alloc = len;
        poly->coeffs = (fq_nmod_struct *) flint_realloc(poly->coeffs,
                                             new_alloc*sizeof(fq_nmod_struct));
        for (i = poly->coeff_alloc; i < new_alloc; i++)
        {
            fq_nmod_init(poly->coeffs + i, dctx->fqctx);
        }        
        poly->coeff_alloc = new_alloc;
    }
}

void fq_nmod_mpolyd_set_nvars(fq_nmod_mpolyd_t poly, slong nvars) {

    poly->nvars = nvars;
    if (poly->degb_alloc < nvars) {
        poly->deg_bounds = (slong *) flint_realloc(poly->deg_bounds, nvars
                                                               *sizeof(slong));
        poly->degb_alloc = nvars;
    }
}

void fq_nmod_mpolyd_zero(fq_nmod_mpolyd_t poly, const fq_nmod_mpolyd_ctx_t dctx)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    fq_nmod_set_ui(poly->coeffs + 0, UWORD(0), dctx->fqctx);
}

void fq_nmod_mpolyd_clear(fq_nmod_mpolyd_t poly, const fq_nmod_mpolyd_ctx_t dctx)
{
    slong i;
    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fq_nmod_clear(poly->coeffs + i, dctx->fqctx);
    }
    flint_free(poly->deg_bounds);
    flint_free(poly->coeffs);
    poly->deg_bounds = NULL;
    poly->coeffs = NULL;
}

void nmod_mpoly_convert_to_fq_nmod_mpolyd(
                       fq_nmod_mpolyd_t poly1, const fq_nmod_mpolyd_ctx_t dctx,
                          const nmod_mpoly_t poly2, const nmod_mpoly_ctx_t ctx)
{
    slong degb_prod;
    slong i, j, N;
    slong * exps;
    const slong * perm = dctx->perm;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    fq_nmod_mpolyd_set_nvars(poly1, ctx->minfo->nvars);

    FLINT_ASSERT(poly2->bits <= FLINT_BITS);

    if (poly2->length == 0)
    {
        fq_nmod_mpolyd_zero(poly1, dctx);
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

    fq_nmod_mpolyd_fit_length(poly1, degb_prod, dctx);
    for (i = 0; i < degb_prod; i++)
    {
        fq_nmod_zero(poly1->coeffs + i, dctx->fqctx);
    }

    N = mpoly_words_per_exp(poly2->bits, ctx->minfo);
    for (i = 0; i < poly2->length; i++)
    {
        slong off = 0;

        mpoly_get_monomial_ui((ulong *)exps, poly2->exps + N*i,
                                                      poly2->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            off = exps[perm[j]] + poly1->deg_bounds[j]*off;
        }
        fq_nmod_set_ui(poly1->coeffs + off, poly2->coeffs[i], dctx->fqctx);
    }

    TMP_END;
}


void nmod_mpoly_convert_from_fq_nmod_mpolyd(
                               nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx,
                     const fq_nmod_mpolyd_t B, const fq_nmod_mpolyd_ctx_t dctx)
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

        if (fq_nmod_is_zero(B->coeffs + i, dctx->fqctx))
            continue;

        for (j = B->nvars - 1; j >= 0; j--) 
        {
            ulong m = B->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            exps[perm[j]] = e;
        }
        FLINT_ASSERT(k == 0);

        /* need special function to convert F_q element to F_p element */
        for (j=1; j < (B->coeffs + i)->length; j++) {
            FLINT_ASSERT((B->coeffs + i)->coeffs[j] == 0);
        }
        nmod_mpoly_set_term_ui_ui(A, (B->coeffs + i)->coeffs[0], exps, ctx);
    }

    TMP_END;
}

void fq_nmod_mpolyd_print(fq_nmod_mpolyd_t poly, const fq_nmod_mpolyd_ctx_t dctx)
{
    int first = 0;
    slong i, j;
    slong degb_prod;

    flint_printf("[ ");
    degb_prod = WORD(1);
    for (j = 0; j < poly->nvars; j++)
    {
        flint_printf("%wd ", poly->deg_bounds[j]);
        degb_prod *= poly->deg_bounds[j];
    }
    flint_printf("]: ");

    first = 1;
    for (i = 0; i < degb_prod; i++)
    {
        ulong k = i;

        if (fq_nmod_is_zero(poly->coeffs + i, dctx->fqctx))
            continue;

        if (!first)
            printf(" + ");

        flint_printf("(");
        fq_nmod_print_pretty(poly->coeffs + i, dctx->fqctx);
        flint_printf(")");

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
    {
        flint_printf("0");
    }
}


void fq_nmod_mpolyd_gcd_brown_univar(fq_nmod_mpolyd_t G,
                             fq_nmod_mpolyd_t Abar,      fq_nmod_mpolyd_t Bbar,
                       const fq_nmod_mpolyd_t A,   const fq_nmod_mpolyd_t B,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{
    slong Alen, Blen;

    Alen = A->deg_bounds[0];
    while (fq_nmod_is_zero(A->coeffs + Alen-1, dctx->fqctx))
    {
        Alen--;
        if (Alen == 0)
            break;
    }

    Blen = B->deg_bounds[0];
    while (fq_nmod_is_zero(B->coeffs + Blen-1, dctx->fqctx))
    {
        Blen--;
        if (Blen == 0)
            break;
    }

    fq_nmod_mpolyd_set_nvars(G, 1);
    fq_nmod_mpolyd_set_nvars(Abar, 1);
    fq_nmod_mpolyd_set_nvars(Bbar, 1);
    if (Alen == 0)
    {
        if (Blen == 0)
        {
            fq_nmod_zero(G->coeffs + 0, dctx->fqctx);
            G->deg_bounds[0] = 1;
            fq_nmod_zero(Abar->coeffs + 0, dctx->fqctx);
            Abar->deg_bounds[0] = 1;
            fq_nmod_zero(Bbar->coeffs + 0, dctx->fqctx);
            Bbar->deg_bounds[0] = 1;
        } else
        {
            fq_nmod_t inv;

            fq_nmod_mpolyd_fit_length(G, Blen, dctx);

            fq_nmod_init(inv, dctx->fqctx);
            fq_nmod_inv(inv, B->coeffs + Blen-1, dctx->fqctx);
            _fq_nmod_poly_scalar_mul_fq_nmod(G->coeffs, B->coeffs, Blen,
                                                         inv, dctx->fqctx);
            fq_nmod_clear(inv, dctx->fqctx);

            G->deg_bounds[0] = Blen;
            fq_nmod_set_ui(Abar->coeffs + 0, 0, dctx->fqctx);
            Abar->deg_bounds[0] = 1;
            fq_nmod_set_ui(Bbar->coeffs + 0, 1, dctx->fqctx);
            Bbar->deg_bounds[0] = 1;
        }

    } else {
        if (Blen == 0)
        {
            fq_nmod_t inv;

            fq_nmod_mpolyd_fit_length(G, Alen, dctx);

            fq_nmod_init(inv, dctx->fqctx);
            fq_nmod_inv(inv, A->coeffs + Alen-1, dctx->fqctx);
            _fq_nmod_poly_scalar_mul_fq_nmod(G->coeffs, A->coeffs, Alen,
                                                         inv, dctx->fqctx);
            fq_nmod_clear(inv, dctx->fqctx);

            G->deg_bounds[0] = Blen;
            fq_nmod_set_ui(Abar->coeffs + 0, 0, dctx->fqctx);
            Abar->deg_bounds[0] = 1;
            fq_nmod_set_ui(Bbar->coeffs + 0, 1, dctx->fqctx);
            Bbar->deg_bounds[0] = 1;

        } else {

            if (Alen >= Blen)
            {
                fq_nmod_t invB;
                fq_nmod_mpolyd_fit_length(G, Blen, dctx);
                fq_nmod_init(invB, dctx->fqctx);
                fq_nmod_inv(invB, B->coeffs + Blen-1, dctx->fqctx);
                G->deg_bounds[0] = _fq_nmod_poly_gcd_euclidean(G->coeffs,
                      A->coeffs, Alen, B->coeffs, Blen, invB, dctx->fqctx);
                fq_nmod_clear(invB, dctx->fqctx);
            } else
            {

                fq_nmod_t invA;
                fq_nmod_mpolyd_fit_length(G, Alen, dctx);
                fq_nmod_init(invA, dctx->fqctx);
                fq_nmod_inv(invA, A->coeffs + Alen-1, dctx->fqctx);
                G->deg_bounds[0] = _fq_nmod_poly_gcd_euclidean(G->coeffs,
                      B->coeffs, Blen, A->coeffs, Alen, invA, dctx->fqctx);
                fq_nmod_clear(invA, dctx->fqctx);
            }

            if (G->deg_bounds[0] <= 1)
            {
                fq_nmod_set_ui(G->coeffs, 1, dctx->fqctx);
            } else
            {
                fq_nmod_t inv;
                fq_nmod_init(inv, dctx->fqctx);
                fq_nmod_inv(inv, G->coeffs + G->deg_bounds[0]-1, dctx->fqctx);
                _fq_nmod_poly_scalar_mul_fq_nmod(G->coeffs, G->coeffs,
                                       G->deg_bounds[0], inv, dctx->fqctx);
                fq_nmod_clear(inv, dctx->fqctx);
            }

            Abar->deg_bounds[0] = Alen - G->deg_bounds[0] + 1;
            fq_nmod_mpolyd_fit_length(Abar, Alen - G->deg_bounds[0] + 1, dctx);

            _fq_nmod_poly_divides(Abar->coeffs, A->coeffs, Alen,
                              G->coeffs, G->deg_bounds[0],
                              G->coeffs + G->deg_bounds[0]-1, dctx->fqctx);

            Bbar->deg_bounds[0] = Blen - G->deg_bounds[0] + 1;
            fq_nmod_mpolyd_fit_length(Bbar, Blen - G->deg_bounds[0] + 1, dctx);
            _fq_nmod_poly_divides(Bbar->coeffs, B->coeffs, Blen,
                              G->coeffs, G->deg_bounds[0],
                              G->coeffs + G->deg_bounds[0]-1, dctx->fqctx);

        }
    }
    return;
}


void fq_nmod_mpolyd_last_content(fq_nmod_poly_t cont, const fq_nmod_mpolyd_t A,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{
    fq_nmod_poly_t temp;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
    slong i, j, Plen;
    slong degb_prod, degb_last=0;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }

    fq_nmod_poly_zero(cont, dctx->fqctx);
    fq_nmod_poly_fit_length(cont, degb_last, dctx->fqctx);

    fq_nmod_poly_init(temp, dctx->fqctx);
    fq_nmod_poly_fit_length(temp, degb_last, dctx->fqctx);

    for (i = 0; i < degb_prod; i += degb_last)
    {
        fq_nmod_struct * P = A->coeffs + i;
        Plen = degb_last;
        while (fq_nmod_is_zero(P + Plen - 1, dctx->fqctx))
        {
            Plen --;
            if (Plen == 0)
                break;
        }

        if (Plen != 0)
        {
            if (cont->length == 0)
            {
                _fq_nmod_poly_make_monic(cont->coeffs, P, Plen, dctx->fqctx);
                cont->length = Plen;
            } else
            {
                fq_nmod_t inv;

                fq_nmod_init(inv, dctx->fqctx);
                if (cont->length < Plen) {
                    fq_nmod_inv(inv, cont->coeffs + cont->length - 1, dctx->fqctx);
                    cont->length = _fq_nmod_poly_gcd_euclidean(temp->coeffs,
                                          P, Plen, cont->coeffs, cont->length,
                                                             inv, dctx->fqctx);
                } else {
                    fq_nmod_inv(inv, P + Plen - 1, dctx->fqctx);
                    cont->length = _fq_nmod_poly_gcd_euclidean(temp->coeffs,
                                          cont->coeffs, cont->length, P, Plen,
                                                             inv, dctx->fqctx);
                }
                fq_nmod_clear(inv, dctx->fqctx);

                if (cont->length == 1)
                {
                    fq_nmod_set_ui(cont->coeffs + 0, UWORD(1), dctx->fqctx);
                } else
                {
                    _fq_nmod_poly_make_monic(cont->coeffs, temp->coeffs,
                                                    cont->length, dctx->fqctx);
                }
            }
        }
    }

    fq_nmod_poly_clear(temp, dctx->fqctx);
}


slong fq_nmod_mpolyd_last_degree(const fq_nmod_mpolyd_t A,
                                               const fq_nmod_mpolyd_ctx_t dctx) 
{                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    slong i, j, Plen, degree;
    slong degb_prod, degb_last=0;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }

    degree = -WORD(1);
    for (i = 0; i < degb_prod; i += degb_last)
    {
        fq_nmod_struct * P = A->coeffs + i;
        Plen = degb_last;
        while (fq_nmod_is_zero(P + Plen - 1, dctx->fqctx))
        {
            Plen --;
            if (Plen == 0)
                break;
        }
        degree = FLINT_MAX(degree, Plen - 1);
    }
    return degree;
}


void fq_nmod_mpolyd_div_last_poly(fq_nmod_mpolyd_t A, fq_nmod_poly_t b,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{
    slong i, j, k;
    slong degb_prod, degb_last=0, new_degb_last;
    fq_nmod_struct * temp;
    TMP_INIT;

    FLINT_ASSERT(b->length != 0);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }
    new_degb_last = degb_last - (b->length - 1);

    TMP_START;
    temp = (fq_nmod_struct *) TMP_ALLOC(new_degb_last*sizeof(fq_nmod_struct));
    for (i = 0; i < new_degb_last; i++)
    {
        fq_nmod_init(temp + i, dctx->fqctx);
    }

    k = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        fq_nmod_struct * P = A->coeffs + i;
        slong Plen = degb_last;
        while (fq_nmod_is_zero(P + Plen - 1, dctx->fqctx))
        {
            Plen --;
            if (Plen == 0)
                break;
        }

        j = 0;
        if (Plen != 0)
        {
            fq_nmod_t inv;

            fq_nmod_init(inv, dctx->fqctx);
            fq_nmod_inv(inv, b->coeffs + b->length - 1, dctx->fqctx);

            FLINT_ASSERT(Plen > b->length - 1);
            _fq_nmod_poly_divides(temp, P, Plen, b->coeffs, b->length, inv, dctx->fqctx);
            while (j < Plen - (b->length - 1))
            {
                fq_nmod_set(A->coeffs + k, temp + j, dctx->fqctx);
                k++;
                j++;
            }

            fq_nmod_clear(inv, dctx->fqctx);
        }
        while (j < new_degb_last)
        {
            fq_nmod_set_ui(A->coeffs + k, UWORD(0), dctx->fqctx);
            k++;
            j++;
        }
    }

    A->deg_bounds[A->nvars - 1] = new_degb_last;

    for (i = 0; i < new_degb_last; i++)
    {
        fq_nmod_clear(temp + i, dctx->fqctx);
    }
    TMP_END;
}


void fq_nmod_mpolyd_mul_last_poly(fq_nmod_mpolyd_t M, fq_nmod_mpolyd_t A,
                             fq_nmod_poly_t b, const fq_nmod_mpolyd_ctx_t dctx)
{
    slong i, j, k;
    slong degb_prod, degb_prod_last=0, degb_last=0, new_degb_last;

    FLINT_ASSERT(M != A);
    FLINT_ASSERT(b->length != 0);

    fq_nmod_mpolyd_set_nvars(M, A->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        M->deg_bounds[j] = A->deg_bounds[j];
        degb_last = A->deg_bounds[j];
        degb_prod_last = degb_prod;
        degb_prod *= degb_last;
    }
    new_degb_last = degb_last + (b->length - 1);

    M->deg_bounds[M->nvars - 1] = degb_last + (b->length - 1);

    fq_nmod_mpolyd_fit_length(M, degb_prod_last*(degb_last + (b->length - 1)), dctx);

    k = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        fq_nmod_struct * P = A->coeffs + i;
        slong Plen = degb_last;
        while (fq_nmod_is_zero(P + Plen - 1, dctx->fqctx))
        {
            Plen --;
            if (Plen == 0)
                break;
        }

        j = 0;
        if (Plen != 0)
        {
            _fq_nmod_poly_mul_classical(M->coeffs + k, P, Plen,
                                            b->coeffs, b->length, dctx->fqctx);
            j +=  Plen + (b->length - 1);
        }
        while (j < new_degb_last)
        {
            fq_nmod_zero(M->coeffs + k + j, dctx->fqctx);
            j++;
        }
        k += new_degb_last;
    }
}


void fq_nmod_mpolyd_mul_scalar(fq_nmod_mpolyd_t A, fq_nmod_t b,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{
    slong j, degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_prod *= A->deg_bounds[j];
    }

    for (j = 0; j < degb_prod; j++)
    {
        fq_nmod_mul(A->coeffs + j, A->coeffs + j, b, dctx->fqctx);
    }
}


void fq_nmod_mpolyd_last_lc(fq_nmod_poly_t lc, const fq_nmod_mpolyd_t A,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{

    slong i, j, lc_len;
    slong degb_prod, degb_last=0;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }

    fq_nmod_poly_zero(lc, dctx->fqctx);

    for (i = degb_prod - degb_last; i >= 0; i -= degb_last)
    {
        lc_len = degb_last;
        while (lc_len > 0)
        {
            if (!fq_nmod_is_zero(A->coeffs + i + lc_len - 1, dctx->fqctx))
            {
                fq_nmod_poly_fit_length(lc, lc_len, dctx->fqctx);
                for (j = 0; j < lc_len; j++)
                {
                    fq_nmod_set(lc->coeffs + j, A->coeffs + i + j, dctx->fqctx);
                }
                lc->length = lc_len;
                return;
            }
            lc_len--;
        }
    }

    return;
}


void fq_nmod_mpolyd_eval_last(fq_nmod_mpolyd_t E, const fq_nmod_mpolyd_t A,
                              fq_nmod_t alpha, const fq_nmod_mpolyd_ctx_t dctx)
{
    slong i, j, k;
    slong degb_prod, degb_last=0, degb_prod_last=0;

    fq_nmod_mpolyd_set_nvars(E, A->nvars - 1);

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

    fq_nmod_mpolyd_fit_length(E, degb_prod_last, dctx);

    j = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        fq_nmod_t pp, v;
        fq_nmod_init(pp, dctx->fqctx);
        fq_nmod_init(v, dctx->fqctx);
        fq_nmod_zero(v, dctx->fqctx);

        k = degb_last;
        while (--k >= 0)
        {
            fq_nmod_mul(pp, v, alpha, dctx->fqctx);
            fq_nmod_add(v, pp, A->coeffs + i + k, dctx->fqctx);
        }
        fq_nmod_set(E->coeffs + j, v, dctx->fqctx);
        j++;

        fq_nmod_clear(pp, dctx->fqctx);
        fq_nmod_clear(v, dctx->fqctx);
    }

    FLINT_ASSERT(j == degb_prod_last);

    return;
}


slong fq_nmod_mpolyd_leadmon(slong * exps, const fq_nmod_mpolyd_t A,
                                               const fq_nmod_mpolyd_ctx_t dctx)
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
        if (!fq_nmod_is_zero(A->coeffs + i, dctx->fqctx))
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


void fq_nmod_mpolyd_set_ui(fq_nmod_mpolyd_t A, mp_limb_t v,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{
    slong j;

    for (j = 0; j < A->nvars; j++)
        A->deg_bounds[j] = WORD(1);

    fq_nmod_set_ui(A->coeffs + 0, v, dctx->fqctx);
}




void fq_nmod_mpolyd_set(fq_nmod_mpolyd_t A, const fq_nmod_mpolyd_t B,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{
    slong j;
    slong degb_prod;

    fq_nmod_mpolyd_set_nvars(A, B->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < B->nvars; j++)
    {
        A->deg_bounds[j] = B->deg_bounds[j];
        degb_prod *= B->deg_bounds[j];
    }

    fq_nmod_mpolyd_fit_length(A, degb_prod, dctx);

    for (j = 0; j < degb_prod; j++)
    {
        fq_nmod_set(A->coeffs + j, B->coeffs + j, dctx->fqctx);
    }
}

void fq_nmod_mpolyd_startinterp(fq_nmod_mpolyd_t fxn,
                    fq_nmod_mpolyd_t newvalue, const fq_nmod_mpolyd_ctx_t dctx)
{
    slong nvars = newvalue->nvars;

    fq_nmod_mpolyd_set_nvars(newvalue, nvars+1);
    newvalue->deg_bounds[nvars] = 1;

    fq_nmod_mpolyd_set(fxn, newvalue, dctx);
}




int fq_nmod_next(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx)
{
    slong i;
    slong deg = fqctx->modulus->length - 1;

    for (i = 0; i < deg; i++) {
        ulong c = nmod_poly_get_coeff_ui(alpha, i);
        c += UWORD(1);
        if (c < fqctx->mod.n) {
            nmod_poly_set_coeff_ui(alpha, i, c);
            return 1;
        }
        nmod_poly_set_coeff_ui(alpha, i, UWORD(0));
    }

    return 0;
}

void fq_nmod_mpolyd_addinterp(fq_nmod_mpolyd_t F, fq_nmod_mpolyd_t T,
             fq_nmod_mpolyd_t N, fq_nmod_poly_t modulus, fq_nmod_t alpha,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{
    slong i, j, k, degb_prod;
    slong nvars = F->nvars;
    fq_nmod_t Nvalue, pp, v, u;
    int changed = 0;
    int carry, Fok, Nok;
    slong * inds, Find, Nind, Tind;
    TMP_INIT;

    fq_nmod_init(Nvalue, dctx->fqctx);
    fq_nmod_init(pp, dctx->fqctx);
    fq_nmod_init(v, dctx->fqctx);
    fq_nmod_init(u, dctx->fqctx);

    FLINT_ASSERT(N->nvars == nvars - 1);
    FLINT_ASSERT(modulus->length > 0);

    fq_nmod_mpolyd_set_nvars(T, nvars);

    degb_prod = 1;
    for (j = 0; j < nvars - 1; j++)
    {
        T->deg_bounds[j] = FLINT_MAX(F->deg_bounds[j], N->deg_bounds[j]);
        degb_prod *= T->deg_bounds[j];
    }
        T->deg_bounds[j] = FLINT_MAX(1 + fq_nmod_mpolyd_last_degree(F, dctx),
                                                  modulus->length);

    fq_nmod_mpolyd_fit_length(T, degb_prod*T->deg_bounds[nvars-1], dctx);

    TMP_START;
    inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (j = 0; j < nvars - 1; j++)
        inds[j] = 0;

    Fok = 1;
    Nok = 1;
    Find = 0;
    Nind = 0;
    Tind = 0;
    for (i = 0; i < degb_prod; i++)
    {
        for (k = 0; k < T->deg_bounds[nvars-1]; k++)
        {
            fq_nmod_zero(T->coeffs + Tind+k, dctx->fqctx);
        }

        fq_nmod_zero(Nvalue, dctx->fqctx);
        if (Nok)
        {
            fq_nmod_set(Nvalue, N->coeffs + Nind, dctx->fqctx);
        }

        if (Fok)
        {
            for (k = 0; k < T->deg_bounds[nvars-1]
                       && k < F->deg_bounds[nvars-1]; k++)
            {
                fq_nmod_set(T->coeffs + Tind+k, F->coeffs + Find+k, dctx->fqctx);
            }

            fq_nmod_zero(v, dctx->fqctx);
            k = F->deg_bounds[nvars-1];
            while (--k >= 0)
            {
                fq_nmod_mul(pp, v, alpha, dctx->fqctx);
                fq_nmod_add(v, pp, F->coeffs + Find + k, dctx->fqctx);
            }

            fq_nmod_sub(v, Nvalue, v, dctx->fqctx);

            if (!fq_nmod_is_zero(v, dctx->fqctx))
            {
                changed = 1;
                for (k = 0; k < modulus->length; k++)
                {
                    fq_nmod_mul(u, modulus->coeffs + k, v, dctx->fqctx);
                    fq_nmod_add(T->coeffs + Tind + k, T->coeffs + Tind + k,
                                                               u, dctx->fqctx);
                }
            }

        } else
        {
            if (!fq_nmod_is_zero(Nvalue, dctx->fqctx))
            {
                changed = 1;
                for (k = 0; k < modulus->length; k++)
                {
                    fq_nmod_mul(T->coeffs + Tind + k, modulus->coeffs + k,
                                                          Nvalue, dctx->fqctx);
                }
            }
        }

        /* move indices to next chunk */

        carry = 1;
        for (j = nvars - 2; j >= 0; j--)
        {
            inds[j] += carry;
            if (inds[j] < T->deg_bounds[j])
            {
                carry = 0;
            } else
            {
                carry = 1;
                inds[j] = 0;
            }
        }

        Tind += T->deg_bounds[nvars-1];

        Find = 0;
        Nind = 0;
        Fok = 1;
        Nok = 1;
        for (j = 0; j < nvars - 1; j++)
        {
            Fok = Fok && (inds[j] < F->deg_bounds[j]);
            Nok = Nok && (inds[j] < N->deg_bounds[j]);
            Find = inds[j] + F->deg_bounds[j]*Find;
            Nind = inds[j] + N->deg_bounds[j]*Nind;
        }
        Find *= F->deg_bounds[nvars - 1];
    }

    if (changed)
    {
        fq_nmod_mpolyd_swap(F, T);
    }

    fq_nmod_clear(Nvalue, dctx->fqctx);
    fq_nmod_clear(pp, dctx->fqctx);
    fq_nmod_clear(v, dctx->fqctx);
    fq_nmod_clear(u, dctx->fqctx);

    TMP_END;
}

int fq_nmod_mpolyd_gcd_brown(fq_nmod_mpolyd_t G,
                               fq_nmod_mpolyd_t Abar, fq_nmod_mpolyd_t Bbar,
                               fq_nmod_mpolyd_t A,    fq_nmod_mpolyd_t B,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{

    int success;
    slong j, bound;
    slong nvars = A->nvars;
    fq_nmod_t alpha, gamma_eval;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, lcA, lcB, gamma;
    fq_nmod_poly_t cGs, cAbars, cBbars, modulus, modulus2;
    fq_nmod_mpolyd_t T, Gs, Abars, Bbars, phiA, phiB, gs, abars, bbars;
    slong leadmon_gs_idx;
    slong * leadmon_gs, * leadmon_Gs;
    slong deggamma, degGs, degA, degB, degAbars, degBbars;

    FLINT_ASSERT(G != A);
    FLINT_ASSERT(G != B);
    FLINT_ASSERT(A->nvars == B->nvars);

    if (A->nvars == 1)
    {
        fq_nmod_mpolyd_gcd_brown_univar(G, Abar, Bbar, A, B, dctx);
        return 1;
    }

    leadmon_gs = (slong *) flint_malloc(nvars*sizeof(slong));
    leadmon_Gs = (slong *) flint_malloc(nvars*sizeof(slong));

    fq_nmod_poly_init(cA, dctx->fqctx);
    fq_nmod_poly_init(cB, dctx->fqctx);
    fq_nmod_mpolyd_last_content(cA, A, dctx);
    fq_nmod_mpolyd_last_content(cB, B, dctx);

    fq_nmod_mpolyd_div_last_poly(A, cA, dctx);
    fq_nmod_mpolyd_div_last_poly(B, cB, dctx);

    fq_nmod_poly_init(cG, dctx->fqctx);
    fq_nmod_poly_gcd_euclidean(cG, cA, cB, dctx->fqctx);

    fq_nmod_poly_init(cAbar, dctx->fqctx);
    fq_nmod_poly_init(cBbar, dctx->fqctx);
    fq_nmod_poly_divides(cAbar, cA, cG, dctx->fqctx);
    fq_nmod_poly_divides(cBbar, cB, cG, dctx->fqctx);

    fq_nmod_poly_init(lcA, dctx->fqctx);
    fq_nmod_poly_init(lcB, dctx->fqctx);
    fq_nmod_mpolyd_last_lc(lcA, A, dctx);
    fq_nmod_mpolyd_last_lc(lcB, B, dctx);

    fq_nmod_poly_init(gamma, dctx->fqctx);
    fq_nmod_poly_gcd_euclidean(gamma, lcA, lcB, dctx->fqctx);

    bound = 1 + fq_nmod_poly_degree(gamma, dctx->fqctx)
               + FLINT_MAX(fq_nmod_mpolyd_last_degree(A, dctx),
                           fq_nmod_mpolyd_last_degree(B, dctx));

    fq_nmod_mpolyd_init(T, nvars, dctx);

    fq_nmod_mpolyd_init(Gs, nvars, dctx);
    fq_nmod_mpolyd_init(Abars, nvars, dctx);
    fq_nmod_mpolyd_init(Bbars, nvars, dctx);
    fq_nmod_mpolyd_init(phiA, nvars - 1, dctx);
    fq_nmod_mpolyd_init(phiB, nvars - 1, dctx);

    fq_nmod_mpolyd_init(gs, nvars - 1, dctx);
    fq_nmod_mpolyd_init(abars, nvars - 1, dctx);
    fq_nmod_mpolyd_init(bbars, nvars - 1, dctx);

    fq_nmod_poly_init(modulus, dctx->fqctx);
    fq_nmod_poly_init(modulus2, dctx->fqctx);

    fq_nmod_init(gamma_eval, dctx->fqctx);
    fq_nmod_init(alpha, dctx->fqctx);

    fq_nmod_poly_one(modulus, dctx->fqctx);

    fq_nmod_set_ui(alpha, 0, dctx->fqctx);
    while (fq_nmod_next(alpha, dctx->fqctx) != 0)
    {
        fq_nmod_poly_evaluate_fq_nmod(gamma_eval, gamma, alpha, dctx->fqctx);

        if (fq_nmod_is_zero(gamma_eval, dctx->fqctx))
            goto break_continue;

        fq_nmod_mpolyd_eval_last(phiA, A, alpha, dctx);
        fq_nmod_mpolyd_eval_last(phiB, B, alpha, dctx);

        success = fq_nmod_mpolyd_gcd_brown(gs, abars, bbars, phiA, phiB, dctx);
        if (success == 0)
            goto break_continue;

        leadmon_gs_idx = fq_nmod_mpolyd_leadmon(leadmon_gs, gs, dctx);

        if (leadmon_gs_idx <= 0)
        {
            FLINT_ASSERT(leadmon_gs_idx == 0);
            fq_nmod_mpolyd_set_ui(Gs, WORD(1), dctx);
            fq_nmod_mpolyd_set(Abars, A, dctx);
            fq_nmod_mpolyd_set(Bbars, B, dctx);
            goto successful;
        }

        if (fq_nmod_poly_degree(modulus, dctx->fqctx) > 0)
        {
            fq_nmod_mpolyd_leadmon(leadmon_Gs, Gs, dctx);
            for (j = 0; j < nvars - 1; j++)
            {
                if (leadmon_gs[j] > leadmon_Gs[j])
                {
                    goto break_continue;
                } else if (leadmon_gs[j] < leadmon_Gs[j])
                {
                    fq_nmod_mpolyd_zero(Gs, dctx);
                    fq_nmod_mpolyd_zero(Abars, dctx);
                    fq_nmod_mpolyd_zero(Bbars, dctx);
                    fq_nmod_poly_one(modulus, dctx->fqctx);
                }
            }
        }

        fq_nmod_mpolyd_mul_scalar(gs, gamma_eval, dctx);

        if (fq_nmod_poly_degree(modulus, dctx->fqctx) > 0)
        {
            fq_nmod_t modulus_eval, inv;

            fq_nmod_init(modulus_eval, dctx->fqctx);
            fq_nmod_init(inv, dctx->fqctx);

            fq_nmod_poly_evaluate_fq_nmod(modulus_eval, modulus, alpha, dctx->fqctx);
            fq_nmod_inv(inv, modulus_eval, dctx->fqctx);

            fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, inv, dctx->fqctx);

            fq_nmod_mpolyd_addinterp(Gs, T, gs, modulus, alpha, dctx);
            fq_nmod_mpolyd_addinterp(Abars, T, abars, modulus, alpha, dctx);
            fq_nmod_mpolyd_addinterp(Bbars, T, bbars, modulus, alpha, dctx);

            fq_nmod_clear(modulus_eval, dctx->fqctx);
            fq_nmod_clear(inv, dctx->fqctx);

        } else
        {
            fq_nmod_poly_one(modulus, dctx->fqctx);
            fq_nmod_mpolyd_startinterp(Gs, gs, dctx);
            fq_nmod_mpolyd_startinterp(Abars, abars, dctx);
            fq_nmod_mpolyd_startinterp(Bbars, bbars, dctx);
        }

        fq_nmod_poly_scalar_mul_fq_nmod(modulus2, modulus, alpha, dctx->fqctx);
        fq_nmod_poly_shift_left(modulus, modulus, 1, dctx->fqctx);
        fq_nmod_poly_sub(modulus, modulus, modulus2, dctx->fqctx);

        if (fq_nmod_poly_degree(modulus, dctx->fqctx) < bound)
            continue;

        deggamma = fq_nmod_poly_degree(gamma, dctx->fqctx);
        degGs = fq_nmod_mpolyd_last_degree(Gs, dctx);
        degA = fq_nmod_mpolyd_last_degree(A, dctx);
        degB = fq_nmod_mpolyd_last_degree(B, dctx);
        degAbars = fq_nmod_mpolyd_last_degree(Abars, dctx);
        degBbars = fq_nmod_mpolyd_last_degree(Bbars, dctx);

        if (   deggamma + degA == degGs + degAbars
            && deggamma + degB == degGs + degBbars
           )
        {
            goto successful;
        } else
        {
            fq_nmod_mpolyd_zero(Gs, dctx);
            fq_nmod_mpolyd_zero(Abars, dctx);
            fq_nmod_mpolyd_zero(Bbars, dctx);
            fq_nmod_poly_one(modulus, dctx->fqctx);
            continue;
        }

break_continue:

        (void)(NULL);
    }

    success = 0;

cleanup:

    flint_free(leadmon_gs);
    flint_free(leadmon_Gs);

    fq_nmod_poly_clear(cA, dctx->fqctx);
    fq_nmod_poly_clear(cB, dctx->fqctx);
    fq_nmod_poly_clear(cG, dctx->fqctx);

    fq_nmod_poly_clear(cAbar, dctx->fqctx);
    fq_nmod_poly_clear(cBbar, dctx->fqctx);

    fq_nmod_poly_clear(lcA, dctx->fqctx);
    fq_nmod_poly_clear(lcB, dctx->fqctx);

    fq_nmod_poly_clear(gamma, dctx->fqctx);

    fq_nmod_mpolyd_clear(T, dctx);

    fq_nmod_mpolyd_clear(Gs, dctx);
    fq_nmod_mpolyd_clear(Abars, dctx);
    fq_nmod_mpolyd_clear(Bbars, dctx);
    fq_nmod_mpolyd_clear(phiA, dctx);
    fq_nmod_mpolyd_clear(phiB, dctx);

    fq_nmod_mpolyd_clear(gs, dctx);
    fq_nmod_mpolyd_clear(abars, dctx);
    fq_nmod_mpolyd_clear(bbars, dctx);

    fq_nmod_poly_clear(modulus, dctx->fqctx);
    fq_nmod_poly_clear(modulus2, dctx->fqctx);

    fq_nmod_clear(gamma_eval, dctx->fqctx);
    fq_nmod_clear(alpha, dctx->fqctx);

    return success;

successful:

    fq_nmod_poly_init(cGs, dctx->fqctx);
    fq_nmod_poly_init(cAbars, dctx->fqctx);
    fq_nmod_poly_init(cBbars, dctx->fqctx);
    fq_nmod_mpolyd_last_content(cGs, Gs, dctx);
    fq_nmod_mpolyd_last_content(cAbars, Abars, dctx);
    fq_nmod_mpolyd_last_content(cBbars, Bbars, dctx);

    fq_nmod_mpolyd_div_last_poly(Gs, cGs, dctx);
    fq_nmod_mpolyd_div_last_poly(Abars, cAbars, dctx);
    fq_nmod_mpolyd_div_last_poly(Bbars, cBbars, dctx);

    fq_nmod_mpolyd_mul_last_poly(G, Gs, cG, dctx);
    fq_nmod_mpolyd_mul_last_poly(Abar, Abars, cAbar, dctx);
    fq_nmod_mpolyd_mul_last_poly(Bbar, Bbars, cBbar, dctx);

    fq_nmod_poly_clear(cGs, dctx->fqctx);
    fq_nmod_poly_clear(cAbars, dctx->fqctx);
    fq_nmod_poly_clear(cBbars, dctx->fqctx);

    success = 1;
    goto cleanup;
}


/*
    helper functions for
    large prime version of the modular gcd
*/

void nmod_mpolyd_lgprime_eval_last(fq_nmod_mpolyd_t E, const fq_nmod_mpolyd_ctx_t efctx,
                                     nmod_mpolyd_t A, const nmodf_ctx_t fctx)
{
    slong i, j, k;
    slong degb_prod, degb_last=0, degb_prod_last=0;
    nmod_poly_t P;

    fq_nmod_mpolyd_set_nvars(E, A->nvars - 1);

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

    fq_nmod_mpolyd_fit_length(E, degb_prod_last, efctx);

    P->mod.n = fctx->mod.n;
    P->mod.ninv = fctx->mod.ninv;
    P->mod.norm = fctx->mod.norm;
    P->alloc = degb_last;

    j = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        k = degb_last;
        while (k > 0 && (A->coeffs + i)[k-1] == UWORD(0))
            k--;
        P->coeffs = A->coeffs + i;
        P->length = k;
        nmod_poly_rem(E->coeffs + j, P, efctx->fqctx->modulus);
        j++;
    }

    FLINT_ASSERT(j == degb_prod_last);

    return;
}

void nmod_mpolyd_lgprime_startinterp(nmod_mpolyd_t A, const nmodf_ctx_t fctx,
                         fq_nmod_mpolyd_t E, const fq_nmod_mpolyd_ctx_t efctx)
{

    slong i, j, k;
    slong degb_prod, degb_last;

    nmod_mpolyd_set_nvars(A, E->nvars + 1);

    degb_prod = WORD(1);
    for (j = 0; j < E->nvars; j++)
    {
        A->deg_bounds[j] = E->deg_bounds[j];
        degb_prod *= E->deg_bounds[j];
    }
    degb_last = nmod_poly_degree(efctx->fqctx->modulus);
    A->deg_bounds[E->nvars] = degb_last;
    degb_prod *= degb_last;

    nmod_mpolyd_fit_length(A, degb_prod);

    j = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        FLINT_ASSERT((E->coeffs + j)->length <= degb_last);
        for (k = 0; k < (E->coeffs + j)->length; k++)
            (A->coeffs + i)[k] = (E->coeffs + j)->coeffs[k];

        for (; k < degb_last; k++)
            (A->coeffs + i)[k] = UWORD(0);

        j++;
    }

    return;
}


void nmod_mpolyd_lgprime_addinterp(nmod_mpolyd_t F, nmod_mpolyd_t T,
                              const nmodf_ctx_t fctx, fq_nmod_mpolyd_t N,
                                 nmod_poly_t modulus, fq_nmod_t modulus_eval,
                                              const fq_nmod_mpolyd_ctx_t efctx)
{


    slong i, j, k, degb_prod;
    slong nvars = F->nvars;
    fq_nmod_t Nvalue, v, u;
    int changed = 0;
    int carry, Fok, Nok;
    slong * inds, Find, Nind, Tind;
    nmod_poly_t P;
    TMP_INIT;

    fq_nmod_init(Nvalue, efctx->fqctx);
    fq_nmod_init(v, efctx->fqctx);
    fq_nmod_init(u, efctx->fqctx);

    FLINT_ASSERT(N->nvars == nvars - 1);
    FLINT_ASSERT(modulus->length > 0);

    nmod_mpolyd_set_nvars(T, nvars);

    degb_prod = 1;
    for (j = 0; j < nvars - 1; j++)
    {
        T->deg_bounds[j] = FLINT_MAX(F->deg_bounds[j], N->deg_bounds[j]);
        degb_prod *= T->deg_bounds[j];
    }
        T->deg_bounds[j] = FLINT_MAX(1 + nmod_mpolyd_last_degree(F, fctx),
                                       modulus->length
                                       + nmod_poly_degree(efctx->fqctx->modulus)
                                    );
    
    nmod_mpolyd_fit_length(T, degb_prod*T->deg_bounds[nvars-1]);

    TMP_START;
    inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (j = 0; j < nvars - 1; j++)
        inds[j] = 0;

    P->mod.n = fctx->mod.n;
    P->mod.ninv = fctx->mod.ninv;
    P->mod.norm = fctx->mod.norm;
    P->alloc = F->deg_bounds[nvars-1];

    Fok = 1;
    Nok = 1;
    Find = 0;
    Nind = 0;
    Tind = 0;
    for (i = 0; i < degb_prod; i++)
    {
        for (k = 0; k < T->deg_bounds[nvars-1]; k++)
        {
            T->coeffs[Tind+k] = UWORD(0);
        }

        fq_nmod_zero(Nvalue, efctx->fqctx);
        if (Nok)
        {
            fq_nmod_set(Nvalue, N->coeffs + Nind, efctx->fqctx);
        }

        if (Fok)
        {
            for (k = 0; k < T->deg_bounds[nvars-1]
                                         && k < F->deg_bounds[nvars-1]; k++)
            {
                T->coeffs[Tind+k] = F->coeffs[Find+k];
            }
            k = F->deg_bounds[nvars-1];
            while (k > 0 && (F->coeffs[Find+k-1] == UWORD(0)))
                k--;
            P->coeffs = F->coeffs + Find;
            P->length = k;
            nmod_poly_rem(v, P, efctx->fqctx->modulus);
            fq_nmod_sub(Nvalue, Nvalue, v, efctx->fqctx);
        }

        if (!fq_nmod_is_zero(Nvalue, efctx->fqctx))
        {
            changed = 1;

            fq_nmod_mul(u, Nvalue, modulus_eval, efctx->fqctx);
            nmod_poly_mul(v, u, modulus);
            FLINT_ASSERT(v->length <= T->deg_bounds[nvars-1]);
            for (k = 0; k < v->length; k++)
                T->coeffs[Tind + k] = nmod_add(T->coeffs[Tind + k],
                                                      v->coeffs[k], fctx->mod);
        }

        /* move indices to next chunk */

        carry = 1;
        for (j = nvars - 2; j >= 0; j--)
        {
            inds[j] += carry;
            if (inds[j] < T->deg_bounds[j])
            {
                carry = 0;
            } else
            {
                carry = 1;
                inds[j] = 0;
            }
        }

        Tind += T->deg_bounds[nvars-1];

        Find = 0;
        Nind = 0;
        Fok = 1;
        Nok = 1;
        for (j = 0; j < nvars - 1; j++)
        {
            Fok = Fok && (inds[j] < F->deg_bounds[j]);
            Nok = Nok && (inds[j] < N->deg_bounds[j]);
            Find = inds[j] + F->deg_bounds[j]*Find;
            Nind = inds[j] + N->deg_bounds[j]*Nind;
        }
        Find *= F->deg_bounds[nvars - 1];
    }

    if (changed) 
    {
        nmod_mpolyd_swap(F, T);
    }

    fq_nmod_clear(Nvalue, efctx->fqctx);
    fq_nmod_clear(v, efctx->fqctx);
    fq_nmod_clear(u, efctx->fqctx);

    TMP_END;
    return;
}





int nmod_mpolyd_gcd_brown_lgprime(nmod_mpolyd_t G,
                                 nmod_mpolyd_t Abar,  nmod_mpolyd_t Bbar,
                                 nmod_mpolyd_t A, nmod_mpolyd_t B,
                                                       const nmodf_ctx_t fctx)
{

    int success;
    slong j, bound, deg;
    slong nvars = A->nvars;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, lcA, lcB, gamma;
    nmod_poly_t cGs, cAbars, cBbars, modulus, modulus_eval;
    nmod_mpolyd_t T, Gs, Abars, Bbars;
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

    nmod_mpolyd_init(T, nvars);

    nmod_mpolyd_init(Gs, nvars);
    nmod_mpolyd_init(Abars, nvars);
    nmod_mpolyd_init(Bbars, nvars);

    nmod_poly_init(modulus_eval, fctx->mod.n);
    nmod_poly_init(modulus, fctx->mod.n);
    nmod_poly_one(modulus);

    /* find a starting degree for the ffield extension */
    deg = 0;
    for (j = 0; j<nvars; j++)
    {
        deg +=   FLINT_BIT_COUNT(A->deg_bounds[j])
               + FLINT_BIT_COUNT(B->deg_bounds[j]);
    }
    deg = 2 + deg/FLINT_BIT_COUNT(fctx->mod.n);

    for (; deg < 100; deg++)
    {
        fq_nmod_mpolyd_ctx_t efctx;
        fq_nmod_mpolyd_t phiA, phiB, gs, abars, bbars;
        fq_nmod_t gamma_eval;

        fq_nmod_mpolyd_ctx_init(efctx, nvars - 1, fctx->mod.n, deg);

        fq_nmod_mpolyd_init(phiA, nvars - 1, efctx);
        fq_nmod_mpolyd_init(phiB, nvars - 1, efctx);
        fq_nmod_mpolyd_init(gs, nvars - 1, efctx);
        fq_nmod_mpolyd_init(abars, nvars - 1, efctx);
        fq_nmod_mpolyd_init(bbars, nvars - 1, efctx);
        fq_nmod_init(gamma_eval, efctx->fqctx);

        nmod_poly_rem(gamma_eval, gamma, efctx->fqctx->modulus);
        if (fq_nmod_is_zero(gamma_eval, efctx->fqctx))
            goto break_continue;

        nmod_mpolyd_lgprime_eval_last(phiA, efctx, A, fctx);
        nmod_mpolyd_lgprime_eval_last(phiB, efctx, B, fctx);

        success = fq_nmod_mpolyd_gcd_brown(gs, abars, bbars, phiA, phiB, efctx);
        if (success == 0)
            goto break_continue;

        leadmon_gs_idx = fq_nmod_mpolyd_leadmon(leadmon_gs, gs, efctx);

        if (leadmon_gs_idx <= 0)
        {
            FLINT_ASSERT(leadmon_gs_idx == 0);
            nmod_mpolyd_set_ui(Gs, WORD(1));
            nmod_mpolyd_set(Abars, A);
            nmod_mpolyd_set(Bbars, B);
            fq_nmod_mpolyd_clear(phiA, efctx);
            fq_nmod_mpolyd_clear(phiB, efctx);
            fq_nmod_mpolyd_clear(gs, efctx);
            fq_nmod_mpolyd_clear(abars, efctx);
            fq_nmod_mpolyd_clear(bbars, efctx);
            fq_nmod_clear(gamma_eval, efctx->fqctx);
            fq_nmod_mpolyd_ctx_clear(efctx);
            goto successful;
        }

        if (nmod_poly_degree(modulus) > 0)
        {
            nmod_mpolyd_leadmon(leadmon_Gs, Gs);
            for (j = 0; j < nvars - 1; j++)
            {
                if (leadmon_gs[j] > leadmon_Gs[j])
                {
                    goto break_continue;

                } else if (leadmon_gs[j] < leadmon_Gs[j])
                {
                    nmod_poly_one(modulus);
                }
            }
        }

        fq_nmod_mpolyd_mul_scalar(gs, gamma_eval, efctx);

        if (nmod_poly_degree(modulus) > 0)
        {
            nmod_poly_rem(modulus_eval, modulus, efctx->fqctx->modulus);
            FLINT_ASSERT(!fq_nmod_is_zero(modulus_eval, efctx->fqctx));
            fq_nmod_inv(modulus_eval, modulus_eval, efctx->fqctx);

            nmod_mpolyd_lgprime_addinterp(Gs,    T, fctx, gs,    modulus, modulus_eval, efctx);
            nmod_mpolyd_lgprime_addinterp(Abars, T, fctx, abars, modulus, modulus_eval, efctx);
            nmod_mpolyd_lgprime_addinterp(Bbars, T, fctx, bbars, modulus, modulus_eval, efctx);

        } else
        {
            nmod_poly_one(modulus);
            nmod_mpolyd_lgprime_startinterp(Gs, fctx, gs, efctx);
            nmod_mpolyd_lgprime_startinterp(Abars, fctx, abars, efctx);
            nmod_mpolyd_lgprime_startinterp(Bbars, fctx, bbars, efctx);
        }

        nmod_poly_mul(modulus, modulus, efctx->fqctx->modulus);

        if (nmod_poly_degree(modulus) < bound)
            goto break_continue;

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
            fq_nmod_mpolyd_clear(phiA, efctx);
            fq_nmod_mpolyd_clear(phiB, efctx);
            fq_nmod_mpolyd_clear(gs, efctx);
            fq_nmod_mpolyd_clear(abars, efctx);
            fq_nmod_mpolyd_clear(bbars, efctx);
            fq_nmod_clear(gamma_eval, efctx->fqctx);
            fq_nmod_mpolyd_ctx_clear(efctx);
            goto successful;

        } else
        {
            nmod_poly_one(modulus);
            goto break_continue;
        }

break_continue:

        fq_nmod_mpolyd_clear(phiA, efctx);
        fq_nmod_mpolyd_clear(phiB, efctx);
        fq_nmod_mpolyd_clear(gs, efctx);
        fq_nmod_mpolyd_clear(abars, efctx);
        fq_nmod_mpolyd_clear(bbars, efctx);
        fq_nmod_clear(gamma_eval, efctx->fqctx);
        fq_nmod_mpolyd_ctx_clear(efctx);
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

    nmod_mpolyd_clear(T);
    nmod_mpolyd_clear(Gs);
    nmod_mpolyd_clear(Abars);
    nmod_mpolyd_clear(Bbars);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus_eval);

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
