/*
    Copyright (C) 2018, 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


/*
    permute the variable ording in the dense represenation
    to something favourable for gcd_brown
*/
int fq_nmod_mpolyd_ctx_set_for_gcd(fq_nmod_mpolyd_ctx_t dctx,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
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
    fq_nmod_mpoly_degrees_si(Aexps, A, ctx);
    fq_nmod_mpoly_degrees_si(Bexps, B, ctx);

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
        }
        else
        {
            fq_nmod_t inv;

            fq_nmod_mpolyd_fit_length(G, Blen, dctx->fqctx);

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
    }
    else
    {
        if (Blen == 0)
        {
            fq_nmod_t inv;

            fq_nmod_mpolyd_fit_length(G, Alen, dctx->fqctx);

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
        }
        else
        {
            if (Alen >= Blen)
            {
                fq_nmod_t invB;
                fq_nmod_mpolyd_fit_length(G, Blen, dctx->fqctx);
                fq_nmod_init(invB, dctx->fqctx);
                fq_nmod_inv(invB, B->coeffs + Blen-1, dctx->fqctx);
                G->deg_bounds[0] = _fq_nmod_poly_gcd_euclidean(G->coeffs,
                      A->coeffs, Alen, B->coeffs, Blen, invB, dctx->fqctx);
                fq_nmod_clear(invB, dctx->fqctx);
            }
            else
            {

                fq_nmod_t invA;
                fq_nmod_mpolyd_fit_length(G, Alen, dctx->fqctx);
                fq_nmod_init(invA, dctx->fqctx);
                fq_nmod_inv(invA, A->coeffs + Alen-1, dctx->fqctx);
                G->deg_bounds[0] = _fq_nmod_poly_gcd_euclidean(G->coeffs,
                      B->coeffs, Blen, A->coeffs, Alen, invA, dctx->fqctx);
                fq_nmod_clear(invA, dctx->fqctx);
            }

            if (G->deg_bounds[0] <= 1)
            {
                fq_nmod_set_ui(G->coeffs, 1, dctx->fqctx);
            }
            else
            {
                fq_nmod_t inv;
                fq_nmod_init(inv, dctx->fqctx);
                fq_nmod_inv(inv, G->coeffs + G->deg_bounds[0]-1, dctx->fqctx);
                _fq_nmod_poly_scalar_mul_fq_nmod(G->coeffs, G->coeffs,
                                       G->deg_bounds[0], inv, dctx->fqctx);
                fq_nmod_clear(inv, dctx->fqctx);
            }

            Abar->deg_bounds[0] = Alen - G->deg_bounds[0] + 1;
            fq_nmod_mpolyd_fit_length(Abar, Alen - G->deg_bounds[0] + 1, dctx->fqctx);

            _fq_nmod_poly_divides(Abar->coeffs, A->coeffs, Alen,
                              G->coeffs, G->deg_bounds[0],
                              G->coeffs + G->deg_bounds[0]-1, dctx->fqctx);

            Bbar->deg_bounds[0] = Blen - G->deg_bounds[0] + 1;
            fq_nmod_mpolyd_fit_length(Bbar, Blen - G->deg_bounds[0] + 1, dctx->fqctx);
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

    fq_nmod_mpolyd_fit_length(M, degb_prod_last*(degb_last + (b->length - 1)), dctx->fqctx);

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

    fq_nmod_mpolyd_fit_length(E, degb_prod_last, dctx->fqctx);

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

    fq_nmod_mpolyd_fit_length(A, degb_prod, dctx->fqctx);

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

    fq_nmod_mpolyd_fit_length(T, degb_prod*T->deg_bounds[nvars-1], dctx->fqctx);

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

int fq_nmod_mpolyd_gcd_brown_smprime(fq_nmod_mpolyd_t G,
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

    fq_nmod_mpolyd_init(T, nvars, dctx->fqctx);

    fq_nmod_mpolyd_init(Gs, nvars, dctx->fqctx);
    fq_nmod_mpolyd_init(Abars, nvars, dctx->fqctx);
    fq_nmod_mpolyd_init(Bbars, nvars, dctx->fqctx);
    fq_nmod_mpolyd_init(phiA, nvars - 1, dctx->fqctx);
    fq_nmod_mpolyd_init(phiB, nvars - 1, dctx->fqctx);

    fq_nmod_mpolyd_init(gs, nvars - 1, dctx->fqctx);
    fq_nmod_mpolyd_init(abars, nvars - 1, dctx->fqctx);
    fq_nmod_mpolyd_init(bbars, nvars - 1, dctx->fqctx);

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

        success = fq_nmod_mpolyd_gcd_brown_smprime(gs, abars, bbars, phiA, phiB, dctx);
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
                    fq_nmod_mpolyd_zero(Gs, dctx->fqctx);
                    fq_nmod_mpolyd_zero(Abars, dctx->fqctx);
                    fq_nmod_mpolyd_zero(Bbars, dctx->fqctx);
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
            fq_nmod_mpolyd_zero(Gs, dctx->fqctx);
            fq_nmod_mpolyd_zero(Abars, dctx->fqctx);
            fq_nmod_mpolyd_zero(Bbars, dctx->fqctx);
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

    fq_nmod_mpolyd_clear(T, dctx->fqctx);

    fq_nmod_mpolyd_clear(Gs, dctx->fqctx);
    fq_nmod_mpolyd_clear(Abars, dctx->fqctx);
    fq_nmod_mpolyd_clear(Bbars, dctx->fqctx);
    fq_nmod_mpolyd_clear(phiA, dctx->fqctx);
    fq_nmod_mpolyd_clear(phiB, dctx->fqctx);

    fq_nmod_mpolyd_clear(gs, dctx->fqctx);
    fq_nmod_mpolyd_clear(abars, dctx->fqctx);
    fq_nmod_mpolyd_clear(bbars, dctx->fqctx);

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
