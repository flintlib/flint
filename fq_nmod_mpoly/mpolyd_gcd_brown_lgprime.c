/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


void fq_nmod_mpolyd_lgprime_eval_last(
    fq_nmod_mpolyd_t E,
    const fq_nmod_mpolyd_ctx_t edctx,
    fq_nmod_mpolyd_t A,
    const fq_nmod_embed_t emb)
{
    slong i, j, k;
    slong degb_prod, degb_last=0, degb_prod_last=0;
    fq_nmod_poly_t P;

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

    fq_nmod_mpolyd_fit_length(E, degb_prod_last, edctx);

    P->alloc = degb_last;

    j = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        k = degb_last;
        while (k > 0 && fq_nmod_is_zero(A->coeffs + i + k-1, edctx->fqctx))
            k--;
        P->coeffs = A->coeffs + i;
        P->length = k;
        _fq_nmod_embed_sm_to_lg(E->coeffs + j, P, emb);
        j++;
    }

    FLINT_ASSERT(j == degb_prod_last);

    return;
}


void fq_nmod_mpolyd_lgprime_startinterp(
    fq_nmod_mpolyd_t A,
    const fq_nmod_mpolyd_ctx_t dctx,
    fq_nmod_mpolyd_t E,
    const fq_nmod_embed_t emb)
{
    slong i, j, k;
    slong degb_prod, degb_last;
    fq_nmod_poly_t t;

    fq_nmod_mpolyd_set_nvars(A, E->nvars + 1);

    degb_prod = WORD(1);
    for (j = 0; j < E->nvars; j++)
    {
        A->deg_bounds[j] = E->deg_bounds[j];
        degb_prod *= E->deg_bounds[j];
    }
    degb_last = fq_nmod_poly_degree(emb->h, dctx->fqctx);
    A->deg_bounds[E->nvars] = degb_last;
    degb_prod *= degb_last;

    fq_nmod_mpolyd_fit_length(A, degb_prod, dctx);

    fq_nmod_poly_init(t, dctx->fqctx);

    j = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        _fq_nmod_embed_lg_to_sm(t, E->coeffs + j, emb);

        FLINT_ASSERT(t->length <= degb_last);
        for (k = 0; k < t->length; k++)
            fq_nmod_set(A->coeffs + i + k, t->coeffs + k, dctx->fqctx);

        for (; k < degb_last; k++)
            fq_nmod_zero(A->coeffs + i + k, dctx->fqctx);

        j++;
    }

    fq_nmod_poly_clear(t, dctx->fqctx);

    return;
}


void fq_nmod_mpolyd_lgprime_addinterp(
    fq_nmod_mpolyd_t F,
    fq_nmod_mpolyd_t T,
    const fq_nmod_mpolyd_ctx_t dctx,
    fq_nmod_mpolyd_t N,
    const fq_nmod_mpolyd_ctx_t edctx,
    fq_nmod_poly_t modulus,
    fq_nmod_t modulus_eval,
    const fq_nmod_embed_t emb)
{
    slong i, j, k, degb_prod, Tlast_degb;
    slong nvars = F->nvars;
    fq_nmod_t Nvalue, v, u;
    fq_nmod_poly_t u_sm, v_sm;
    int changed = 0;
    int carry, Fok, Nok;
    slong * inds, Find, Nind, Tind;
    fq_nmod_poly_t P;
    TMP_INIT;

    fq_nmod_init(Nvalue, emb->lgctx);
    fq_nmod_init(v, emb->lgctx);
    fq_nmod_init(u, emb->lgctx);
    fq_nmod_poly_init(u_sm, dctx->fqctx);
    fq_nmod_poly_init(v_sm, dctx->fqctx);

    FLINT_ASSERT(N->nvars == nvars - 1);
    FLINT_ASSERT(modulus->length > 0);

    fq_nmod_mpolyd_set_nvars(T, nvars);

    degb_prod = 1;
    for (j = 0; j < nvars - 1; j++)
    {
        T->deg_bounds[j] = FLINT_MAX(F->deg_bounds[j], N->deg_bounds[j]);
        degb_prod *= T->deg_bounds[j];
    }

    Tlast_degb = FLINT_MAX(1 + fq_nmod_mpolyd_last_degree(F, dctx),
                   modulus->length + fq_nmod_poly_degree(emb->h, dctx->fqctx));
    T->deg_bounds[nvars-1] = Tlast_degb;
    fq_nmod_mpolyd_fit_length(T, degb_prod*Tlast_degb, dctx);

    TMP_START;
    inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (j = 0; j < nvars - 1; j++)
        inds[j] = 0;


    P->alloc = F->deg_bounds[nvars-1];

    Fok = 1;
    Nok = 1;
    Find = 0;
    Nind = 0;
    Tind = 0;
    for (i = 0; i < degb_prod; i++)
    {
        for (k = 0; k < Tlast_degb; k++)
        {
            fq_nmod_zero(T->coeffs + Tind+k, dctx->fqctx);
        }

        fq_nmod_zero(Nvalue, emb->lgctx);
        if (Nok)
        {
            fq_nmod_set(Nvalue, N->coeffs + Nind, emb->lgctx);
        }

        if (Fok)
        {
            for (k = 0; k < Tlast_degb && k < F->deg_bounds[nvars-1]; k++)
            {
                fq_nmod_set(T->coeffs + Tind+k, F->coeffs + Find+k, dctx->fqctx);
            }
            k = F->deg_bounds[nvars-1];
            while (k > 0 && fq_nmod_is_zero(F->coeffs + Find+k-1, dctx->fqctx))
                k--;
            P->coeffs = F->coeffs + Find;
            P->length = k;
            _fq_nmod_embed_sm_to_lg(v, P, emb);
            fq_nmod_sub(Nvalue, Nvalue, v, emb->lgctx);
        }

        if (!fq_nmod_is_zero(Nvalue, emb->lgctx))
        {
            changed = 1;

            fq_nmod_mul(u, Nvalue, modulus_eval, emb->lgctx);
            _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
            fq_nmod_poly_mul(v_sm, u_sm, modulus, dctx->fqctx);
            FLINT_ASSERT(v_sm->length <= T->deg_bounds[nvars-1]);
            for (k = 0; k < v_sm->length; k++)
            {
                fq_nmod_add(T->coeffs + Tind + k, T->coeffs + Tind + k,
                                                v_sm->coeffs + k, dctx->fqctx);
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

        Tind += Tlast_degb;

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

    fq_nmod_clear(Nvalue, emb->lgctx);
    fq_nmod_clear(v, emb->lgctx);
    fq_nmod_clear(u, emb->lgctx);
    fq_nmod_poly_clear(u_sm, dctx->fqctx);
    fq_nmod_poly_clear(v_sm, dctx->fqctx);

    TMP_END;
    return;
}


int fq_nmod_mpolyd_gcd_brown_lgprime(fq_nmod_mpolyd_t G,
                                 fq_nmod_mpolyd_t Abar,  fq_nmod_mpolyd_t Bbar,
                                 fq_nmod_mpolyd_t A, fq_nmod_mpolyd_t B,
                                               const fq_nmod_mpolyd_ctx_t dctx)
{
    int success;
    slong k, j, bound;
    slong nvars = A->nvars;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, lcA, lcB, gamma;
    fq_nmod_poly_t cGs, cAbars, cBbars, modulus;
    fq_nmod_mpolyd_t T, Gs, Abars, Bbars;
    slong leadmon_gs_idx;
    slong * leadmon_gs, * leadmon_Gs;
    slong deggamma, degGs, degA, degB, degAbars, degBbars;
    flint_rand_t state;
    fq_nmod_embed_struct * embed;
    mp_limb_t p = dctx->fqctx->modulus->mod.n;
    slong m = nmod_poly_degree(dctx->fqctx->modulus);
    slong n;

    FLINT_ASSERT(G != A);
    FLINT_ASSERT(G != B);
    FLINT_ASSERT(A->nvars == B->nvars);

    if (A->nvars == 1)
    {
        fq_nmod_mpolyd_gcd_brown_univar(G, Abar, Bbar, A, B, dctx);
        return 1;
    }

    flint_randinit(state);

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

    fq_nmod_poly_init(modulus, dctx->fqctx);
    fq_nmod_poly_one(modulus, dctx->fqctx);

    embed = (fq_nmod_embed_struct *) flint_malloc(m*sizeof(fq_nmod_embed_struct));

    /* n is the degree of the extension */
    n = (FLINT_BITS/2)/(m*FLINT_BIT_COUNT(p));
    for (n = FLINT_MAX(n, WORD(2)); n < 1000; n++)
    {
        nmod_poly_t ext_modulus;
        fq_nmod_ctx_t ext_fqctx;
        fq_nmod_mpolyd_ctx_t edctx;
        fq_nmod_mpolyd_t phiA, phiB, gs, abars, bbars;
        fq_nmod_t gamma_eval;
        fq_nmod_t modulus_eval;

        /* init edctx with modulus of degree m*n */
        nmod_poly_init2(ext_modulus, p, m*n + 1);
        nmod_poly_randtest_sparse_irreducible(ext_modulus, state, m*n + 1);
        fq_nmod_ctx_init_modulus(ext_fqctx, ext_modulus, "$");
        fq_nmod_mpolyd_ctx_init2(edctx, nvars - 1, ext_fqctx);
        fq_nmod_ctx_clear(ext_fqctx);
        nmod_poly_clear(ext_modulus);

        _fq_nmod_embed_array_init(embed, edctx->fqctx, dctx->fqctx);

        fq_nmod_mpolyd_init(phiA, nvars - 1, edctx);
        fq_nmod_mpolyd_init(phiB, nvars - 1, edctx);
        fq_nmod_mpolyd_init(gs, nvars - 1, edctx);
        fq_nmod_mpolyd_init(abars, nvars - 1, edctx);
        fq_nmod_mpolyd_init(bbars, nvars - 1, edctx);

        fq_nmod_init(modulus_eval, edctx->fqctx);
        fq_nmod_init(gamma_eval, edctx->fqctx);

        for (k = 0; k < m; k++)
        {
            /* make sure reduction does not kill both lc */
            _fq_nmod_embed_sm_to_lg(gamma_eval, gamma, embed + k);
            if (fq_nmod_is_zero(gamma_eval, edctx->fqctx))
            {
                goto break_continue;
            }

            fq_nmod_mpolyd_lgprime_eval_last(phiA, edctx, A, embed + k);
            fq_nmod_mpolyd_lgprime_eval_last(phiB, edctx, B, embed + k);

            success = fq_nmod_mpolyd_gcd_brown_smprime(gs, abars, bbars, phiA, phiB, edctx);
            if (success == 0)
            {
                goto break_continue;
            }

            leadmon_gs_idx = fq_nmod_mpolyd_leadmon(leadmon_gs, gs, edctx);

            if (leadmon_gs_idx <= 0)
            {
                FLINT_ASSERT(leadmon_gs_idx == 0);
                fq_nmod_mpolyd_set_ui(Gs, 1, dctx);
                fq_nmod_mpolyd_set(Abars, A, dctx);
                fq_nmod_mpolyd_set(Bbars, B, dctx);
                fq_nmod_mpolyd_clear(phiA, edctx);
                fq_nmod_mpolyd_clear(phiB, edctx);
                fq_nmod_mpolyd_clear(gs, edctx);
                fq_nmod_mpolyd_clear(abars, edctx);
                fq_nmod_mpolyd_clear(bbars, edctx);
                fq_nmod_clear(modulus_eval, edctx->fqctx);
                fq_nmod_clear(gamma_eval, edctx->fqctx);
                _fq_nmod_embed_array_clear(embed, m);
                fq_nmod_mpolyd_ctx_clear(edctx);
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
                    }
                    else if (leadmon_gs[j] < leadmon_Gs[j])
                    {
                        fq_nmod_mpolyd_zero(Gs, dctx);
                        fq_nmod_mpolyd_zero(Abars, dctx);
                        fq_nmod_mpolyd_zero(Bbars, dctx);
                        fq_nmod_poly_one(modulus, dctx->fqctx);
                    }
                }
            }

            fq_nmod_mpolyd_mul_scalar(gs, gamma_eval, edctx);

            if (fq_nmod_poly_degree(modulus, dctx->fqctx) > 0)
            {
                _fq_nmod_embed_sm_to_lg(modulus_eval, modulus, embed + k);
                FLINT_ASSERT(!fq_nmod_is_zero(modulus_eval, edctx->fqctx));
                fq_nmod_inv(modulus_eval, modulus_eval, edctx->fqctx);
                fq_nmod_mpolyd_lgprime_addinterp(Gs,    T, dctx, gs,    edctx, modulus, modulus_eval, embed + k);
                fq_nmod_mpolyd_lgprime_addinterp(Abars, T, dctx, abars, edctx, modulus, modulus_eval, embed + k);
                fq_nmod_mpolyd_lgprime_addinterp(Bbars, T, dctx, bbars, edctx, modulus, modulus_eval, embed + k);
            }
            else
            {
                fq_nmod_poly_one(modulus, dctx->fqctx);
                fq_nmod_mpolyd_lgprime_startinterp(Gs, dctx, gs, embed + k);
                fq_nmod_mpolyd_lgprime_startinterp(Abars, dctx, abars, embed + k);
                fq_nmod_mpolyd_lgprime_startinterp(Bbars, dctx, bbars, embed + k);
            }

            fq_nmod_poly_mul(modulus, modulus, (embed + k)->h, dctx->fqctx);

            if (fq_nmod_poly_degree(modulus, dctx->fqctx) < bound)
                goto break_continue;

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
                fq_nmod_mpolyd_clear(phiA, edctx);
                fq_nmod_mpolyd_clear(phiB, edctx);
                fq_nmod_mpolyd_clear(gs, edctx);
                fq_nmod_mpolyd_clear(abars, edctx);
                fq_nmod_mpolyd_clear(bbars, edctx);
                fq_nmod_clear(modulus_eval, edctx->fqctx);
                fq_nmod_clear(gamma_eval, edctx->fqctx);
                _fq_nmod_embed_array_clear(embed, m);
                fq_nmod_mpolyd_ctx_clear(edctx);
                goto successful;
            }
            else
            {
                fq_nmod_poly_one(modulus, dctx->fqctx);
                continue;
            }

    break_continue:
            (void)(NULL);
        }

        fq_nmod_mpolyd_clear(phiA, edctx);
        fq_nmod_mpolyd_clear(phiB, edctx);
        fq_nmod_mpolyd_clear(gs, edctx);
        fq_nmod_mpolyd_clear(abars, edctx);
        fq_nmod_mpolyd_clear(bbars, edctx);
        fq_nmod_clear(modulus_eval, edctx->fqctx);
        fq_nmod_clear(gamma_eval, edctx->fqctx);
        _fq_nmod_embed_array_clear(embed, m);
        fq_nmod_mpolyd_ctx_clear(edctx);
    }

    success = 0;

cleanup:

    flint_free(leadmon_gs);
    flint_free(leadmon_Gs);
    flint_free(embed);

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

    fq_nmod_poly_clear(modulus, dctx->fqctx);

    flint_randclear(state);

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
