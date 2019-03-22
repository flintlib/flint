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

void nmod_mpolyd_lgprime_eval_last(fq_nmod_mpolyd_t E, const fq_nmod_mpolyd_ctx_t efctx,
                                     nmod_mpolyd_t A, const nmod_ctx_t fctx)
{
    slong i, j, k;
    slong degb_prod, degb_last=0, degb_prod_last=0;
    nmod_polydr_t P;

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

    fq_nmod_mpolyd_fit_length(E, degb_prod_last, efctx->fqctx);

    P->alloc = degb_last;

    j = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        k = degb_last;
        while (k > 0 && (A->coeffs + i)[k-1] == UWORD(0))
            k--;
        P->coeffs = A->coeffs + i;
        P->length = k;
        nmod_polydr_rem(E->coeffs + j, P, efctx->fqctx->modulus, fctx);
        j++;
    }

    FLINT_ASSERT(j == degb_prod_last);

    return;
}

void nmod_mpolyd_lgprime_startinterp(nmod_mpolyd_t A, const nmod_ctx_t fctx,
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
    degb_last = nmod_polydr_degree(efctx->fqctx->modulus, fctx);
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
                              const nmod_ctx_t fctx, fq_nmod_mpolyd_t N,
                                 nmod_polydr_t modulus, fq_nmod_t modulus_eval,
                                              const fq_nmod_mpolyd_ctx_t efctx)
{


    slong i, j, k, degb_prod;
    slong nvars = F->nvars;
    fq_nmod_t Nvalue, v, u;
    int changed = 0;
    int carry, Fok, Nok;
    slong * inds, Find, Nind, Tind;
    nmod_polydr_t P;
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
                                       + nmod_polydr_degree(efctx->fqctx->modulus, fctx)
                                    );
    
    nmod_mpolyd_fit_length(T, degb_prod*T->deg_bounds[nvars-1]);

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
            nmod_polydr_rem(v, P, efctx->fqctx->modulus, fctx);
            fq_nmod_sub(Nvalue, Nvalue, v, efctx->fqctx);
        }

        if (!fq_nmod_is_zero(Nvalue, efctx->fqctx))
        {
            changed = 1;

            fq_nmod_mul(u, Nvalue, modulus_eval, efctx->fqctx);
            nmod_polydr_mul(v, u, modulus, fctx);
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
                                                       const nmod_ctx_t fctx)
{

    int success;
    slong j, bound, deg;
    slong nvars = A->nvars;
    nmod_polydr_t cA, cB, cG, cAbar, cBbar, lcA, lcB, gamma;
    nmod_polydr_t cGs, cAbars, cBbars, modulus, modulus_eval;
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

    nmod_polydr_init(cA, fctx);
    nmod_polydr_init(cB, fctx);
    nmod_mpolyd_last_content(cA, A, fctx);
    nmod_mpolyd_last_content(cB, B, fctx);

    nmod_mpolyd_div_last_poly(A, cA, fctx);
    nmod_mpolyd_div_last_poly(B, cB, fctx);

    nmod_polydr_init(cG, fctx);
    nmod_polydr_gcd_euclidean(cG, cA, cB, fctx);

    nmod_polydr_init(cAbar, fctx);
    nmod_polydr_init(cBbar, fctx);
    nmod_polydr_div(cAbar, cA, cG, fctx);
    nmod_polydr_div(cBbar, cB, cG, fctx);

    nmod_polydr_init(lcA, fctx);
    nmod_polydr_init(lcB, fctx);
    nmod_mpolyd_last_lc(lcA, A, fctx);
    nmod_mpolyd_last_lc(lcB, B, fctx);

    nmod_polydr_init(gamma, fctx);
    nmod_polydr_gcd_euclidean(gamma, lcA, lcB, fctx);

    bound = 1 + nmod_polydr_degree(gamma, fctx)
               + FLINT_MAX(nmod_mpolyd_last_degree(A, fctx),
                           nmod_mpolyd_last_degree(B, fctx));

    nmod_mpolyd_init(T, nvars);

    nmod_mpolyd_init(Gs, nvars);
    nmod_mpolyd_init(Abars, nvars);
    nmod_mpolyd_init(Bbars, nvars);

    nmod_polydr_init(modulus_eval, fctx);
    nmod_polydr_init(modulus, fctx);
    nmod_polydr_one(modulus, fctx);

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

        fq_nmod_mpolyd_init(phiA, nvars - 1, efctx->fqctx);
        fq_nmod_mpolyd_init(phiB, nvars - 1, efctx->fqctx);
        fq_nmod_mpolyd_init(gs, nvars - 1, efctx->fqctx);
        fq_nmod_mpolyd_init(abars, nvars - 1, efctx->fqctx);
        fq_nmod_mpolyd_init(bbars, nvars - 1, efctx->fqctx);
        fq_nmod_init(gamma_eval, efctx->fqctx);

        nmod_polydr_rem(gamma_eval, gamma, efctx->fqctx->modulus, fctx);
        if (fq_nmod_is_zero(gamma_eval, efctx->fqctx))
            goto break_continue;

        nmod_mpolyd_lgprime_eval_last(phiA, efctx, A, fctx);
        nmod_mpolyd_lgprime_eval_last(phiB, efctx, B, fctx);

        success = fq_nmod_mpolyd_gcd_brown_smprime(gs, abars, bbars, phiA, phiB, efctx);
        if (success == 0)
            goto break_continue;

        leadmon_gs_idx = fq_nmod_mpolyd_leadmon(leadmon_gs, gs, efctx);

        if (leadmon_gs_idx <= 0)
        {
            FLINT_ASSERT(leadmon_gs_idx == 0);
            nmod_mpolyd_set_ui(Gs, WORD(1));
            nmod_mpolyd_set(Abars, A);
            nmod_mpolyd_set(Bbars, B);
            fq_nmod_mpolyd_clear(phiA, efctx->fqctx);
            fq_nmod_mpolyd_clear(phiB, efctx->fqctx);
            fq_nmod_mpolyd_clear(gs, efctx->fqctx);
            fq_nmod_mpolyd_clear(abars, efctx->fqctx);
            fq_nmod_mpolyd_clear(bbars, efctx->fqctx);
            fq_nmod_clear(gamma_eval, efctx->fqctx);
            fq_nmod_mpolyd_ctx_clear(efctx);
            goto successful;
        }

        if (nmod_polydr_degree(modulus, fctx) > 0)
        {
            nmod_mpolyd_leadmon(leadmon_Gs, Gs);
            for (j = 0; j < nvars - 1; j++)
            {
                if (leadmon_gs[j] > leadmon_Gs[j])
                {
                    goto break_continue;

                } else if (leadmon_gs[j] < leadmon_Gs[j])
                {
                    nmod_polydr_one(modulus, fctx);
                }
            }
        }

        fq_nmod_mpolyd_mul_scalar(gs, gamma_eval, efctx);

        if (nmod_polydr_degree(modulus, fctx) > 0)
        {
            nmod_polydr_rem(modulus_eval, modulus, efctx->fqctx->modulus, fctx);
            FLINT_ASSERT(!fq_nmod_is_zero(modulus_eval, efctx->fqctx));
            fq_nmod_inv(modulus_eval, modulus_eval, efctx->fqctx);

            nmod_mpolyd_lgprime_addinterp(Gs,    T, fctx, gs,    modulus, modulus_eval, efctx);
            nmod_mpolyd_lgprime_addinterp(Abars, T, fctx, abars, modulus, modulus_eval, efctx);
            nmod_mpolyd_lgprime_addinterp(Bbars, T, fctx, bbars, modulus, modulus_eval, efctx);
        }
        else
        {
            nmod_polydr_one(modulus, fctx);
            nmod_mpolyd_lgprime_startinterp(Gs, fctx, gs, efctx);
            nmod_mpolyd_lgprime_startinterp(Abars, fctx, abars, efctx);
            nmod_mpolyd_lgprime_startinterp(Bbars, fctx, bbars, efctx);
        }

        nmod_polydr_mul(modulus, modulus, efctx->fqctx->modulus, fctx);

        if (nmod_polydr_degree(modulus, fctx) < bound)
            goto break_continue;

        deggamma = nmod_polydr_degree(gamma, fctx);
        degGs = nmod_mpolyd_last_degree(Gs, fctx);
        degA = nmod_mpolyd_last_degree(A, fctx);
        degB = nmod_mpolyd_last_degree(B, fctx);
        degAbars = nmod_mpolyd_last_degree(Abars, fctx);
        degBbars = nmod_mpolyd_last_degree(Bbars, fctx);

        if (   deggamma + degA == degGs + degAbars
            && deggamma + degB == degGs + degBbars
           )
        {
            fq_nmod_mpolyd_clear(phiA, efctx->fqctx);
            fq_nmod_mpolyd_clear(phiB, efctx->fqctx);
            fq_nmod_mpolyd_clear(gs, efctx->fqctx);
            fq_nmod_mpolyd_clear(abars, efctx->fqctx);
            fq_nmod_mpolyd_clear(bbars, efctx->fqctx);
            fq_nmod_clear(gamma_eval, efctx->fqctx);
            fq_nmod_mpolyd_ctx_clear(efctx);
            goto successful;
        }
        else
        {
            nmod_polydr_one(modulus, fctx);
            goto break_continue;
        }

break_continue:

        fq_nmod_mpolyd_clear(phiA, efctx->fqctx);
        fq_nmod_mpolyd_clear(phiB, efctx->fqctx);
        fq_nmod_mpolyd_clear(gs, efctx->fqctx);
        fq_nmod_mpolyd_clear(abars, efctx->fqctx);
        fq_nmod_mpolyd_clear(bbars, efctx->fqctx);
        fq_nmod_clear(gamma_eval, efctx->fqctx);
        fq_nmod_mpolyd_ctx_clear(efctx);
    }

    success = 0;

cleanup:

    flint_free(leadmon_gs);
    flint_free(leadmon_Gs);

    nmod_polydr_clear(cA, fctx);
    nmod_polydr_clear(cB, fctx);
    nmod_polydr_clear(cG, fctx);

    nmod_polydr_clear(cAbar, fctx);
    nmod_polydr_clear(cBbar, fctx);

    nmod_polydr_clear(lcA, fctx);
    nmod_polydr_clear(lcB, fctx);

    nmod_polydr_clear(gamma, fctx);

    nmod_mpolyd_clear(T);
    nmod_mpolyd_clear(Gs);
    nmod_mpolyd_clear(Abars);
    nmod_mpolyd_clear(Bbars);

    nmod_polydr_clear(modulus, fctx);
    nmod_polydr_clear(modulus_eval, fctx);

    return success;

successful:

    nmod_polydr_init(cGs, fctx);
    nmod_polydr_init(cAbars, fctx);
    nmod_polydr_init(cBbars, fctx);
    nmod_mpolyd_last_content(cGs, Gs, fctx);
    nmod_mpolyd_last_content(cAbars, Abars, fctx);
    nmod_mpolyd_last_content(cBbars, Bbars, fctx);

    nmod_mpolyd_div_last_poly(Gs, cGs, fctx);
    nmod_mpolyd_div_last_poly(Abars, cAbars, fctx);
    nmod_mpolyd_div_last_poly(Bbars, cBbars, fctx);

    nmod_mpolyd_mul_last_poly(G, Gs, cG, fctx);
    nmod_mpolyd_mul_last_poly(Abar, Abars, cAbar, fctx);
    nmod_mpolyd_mul_last_poly(Bbar, Bbars, cBbar, fctx);

    nmod_polydr_clear(cGs, fctx);
    nmod_polydr_clear(cAbars, fctx);
    nmod_polydr_clear(cBbars, fctx);

    success = 1;
    goto cleanup;
}
