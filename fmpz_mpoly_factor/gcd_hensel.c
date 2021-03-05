/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"

int fmpz_mpolyl_gcd_hensel(
    fmpz_mpoly_t G, slong Gdeg, /* upperbound on deg_X(G) */
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    int success, alpha_bits, gamma_is_one;
    const slong n = ctx->minfo->nvars - 1;
    slong i, k;
    flint_bitcnt_t bits = A->bits;
    fmpz * alphas, * prev_alphas;
    fmpz_t q, mu1, mu2;
    fmpq_t mu;
    fmpz_mpoly_struct * Aevals, * Bevals, * Hevals;
    fmpz_mpoly_struct * H; /* points to A, B, or Hevals + n */
    fmpz_mpoly_struct * Glcs, * Hlcs;
    fmpz_mpoly_struct Hfac[2], Htfac[2];
    slong * Hdegs;
    slong Adegx, Bdegx, gdegx;
    fmpz_mpoly_t t1, t2, g, abar, bbar, hbar;
    flint_rand_t state;

    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    flint_randinit(state);

    Hdegs  = FLINT_ARRAY_ALLOC(n + 1, slong);

    Glcs   = FLINT_ARRAY_ALLOC(3*(n + 1), fmpz_mpoly_struct);
    Hlcs   = Glcs + (n + 1);
    Hevals = Hlcs + (n + 1);
    for (i = 0; i < n + 1; i++)
    {
        fmpz_mpoly_init(Glcs + i, ctx);
        fmpz_mpoly_init(Hlcs + i, ctx);
        fmpz_mpoly_init(Hevals + i, ctx);
    }

    alphas = _fmpz_vec_init(2*n);
    prev_alphas = alphas + n;
    Aevals = FLINT_ARRAY_ALLOC(2*(n + 1), fmpz_mpoly_struct);
    Bevals = Aevals + (n + 1);
    for (i = 0; i < n; i++)
    {
        fmpz_mpoly_init(Aevals + i, ctx);
        fmpz_mpoly_init(Bevals + i, ctx);
    }

    fmpz_init(q);
    fmpq_init(mu);
    fmpz_init(mu1);
    fmpz_init(mu2);

    fmpz_mpoly_init(t1, ctx);
    fmpz_mpoly_init(t2, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(abar, ctx);
    fmpz_mpoly_init(bbar, ctx);
    fmpz_mpoly_init(hbar, ctx);

    fmpz_mpoly_init(Hfac + 0, ctx);
    fmpz_mpoly_init(Hfac + 1, ctx);
    fmpz_mpoly_init(Htfac + 0, ctx);
    fmpz_mpoly_init(Htfac + 1, ctx);

    /* init done */

    alpha_bits = 0;

    /* try all zeros first */
    for (i = 0; i < n; i++)
    {
        fmpz_zero(prev_alphas + i);
        fmpz_zero(alphas + i);
    }

    goto got_alpha;

next_alpha:

    alpha_bits++;
    if (alpha_bits > FLINT_BITS/2)
    {
        success = 0;
        goto cleanup;
    }

    success = 0;

    for (i = 0; i < n; i++)
    {
        ulong l = n_randlimb(state);
        ulong mask = UWORD(1) << alpha_bits;
        if (l & mask)
            fmpz_neg_ui(alphas + i, l & (mask - 1));
        else
            fmpz_set_ui(alphas + i, l & (mask - 1));

        success = success || !fmpz_equal(alphas + i, prev_alphas + i);
    }

    if (!success)
        goto next_alpha;

got_alpha:

    /* ensure deg_X do not drop under evaluation */
    Adegx = fmpz_mpoly_degree_si(A, 0, ctx);
    Bdegx = fmpz_mpoly_degree_si(B, 0, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i, i == n - 1 ? A :
                                       Aevals + i + 1, i + 1, alphas + i, ctx);
        fmpz_mpoly_evaluate_one_fmpz(Bevals + i, i == n - 1 ? B :
                                       Bevals + i + 1, i + 1, alphas + i, ctx);
        if (Adegx != fmpz_mpoly_degree_si(Aevals + i, 0, ctx) ||
            Bdegx != fmpz_mpoly_degree_si(Bevals + i, 0, ctx))
        {
            goto next_alpha;
        }
    }

    /* univariate gcd */
    success = fmpz_mpoly_gcd_cofactors(g, abar, bbar, Aevals + 0, Bevals + 0, ctx) &&
              fmpz_mpoly_gcd(t1, g, abar, ctx) &&
              fmpz_mpoly_gcd(t2, g, bbar, ctx);
    if (!success)
        goto cleanup;

    gdegx = fmpz_mpoly_degree_si(g, 0, ctx);

    if (gdegx == 0)
    {
        /* G is trivial */
        fmpz_mpoly_set(Abar, A, ctx);
        fmpz_mpoly_set(Bbar, B, ctx);
        fmpz_mpoly_one(G, ctx);
        success = 1;
        goto cleanup;
    }
    else if (gdegx > Gdeg)
    {
        goto next_alpha;
    }
    else if (gdegx < Gdeg)
    {
        Gdeg = gdegx;
        for (i = 0; i < n; i++)
            fmpz_set(prev_alphas + i, alphas + i);

        goto next_alpha;
    }

    /* the degbound gdegx (== Gdeg) has at least two witnesses now */

    if (gdegx == Adegx)
    {
        if (fmpz_mpoly_divides(Bbar, B, A, ctx))
        {
            fmpz_mpoly_set(G, A, ctx);
            fmpz_mpoly_one(Abar, ctx);
            success = 1;
            goto cleanup;
        }

        goto next_alpha;
    }
    else if (gdegx == Bdegx)
    {
        if (fmpz_mpoly_divides(Abar, A, B, ctx))
        {
            fmpz_mpoly_set(G, B, ctx);
            fmpz_mpoly_one(Bbar, ctx);
            success = 1;
            goto cleanup;
        }

        goto next_alpha;
    }

    FLINT_ASSERT(0 < gdegx && gdegx < FLINT_MIN(Adegx, Bdegx));

    /* set Hlcs[n], Glcs[n] (gamma), H, and Hevals */
    if (fmpz_mpoly_is_one(t1, ctx))
    {
        fmpz_one(mu1);
        fmpz_zero(mu2);

        fmpz_mpoly_swap(hbar, abar, ctx);

        fmpz_mpolyl_lead_coeff(Hlcs + n, A, 1, ctx);
        fmpz_mpolyl_lead_coeff(t2, B, 1, ctx);
        success = fmpz_mpoly_gcd(Glcs + n, Hlcs + n, t2, ctx);
        if (!success)
            goto cleanup;

        H = (fmpz_mpoly_struct *) A;

        gamma_is_one = fmpz_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                fmpz_mpoly_swap(Hevals + i, Aevals + i, ctx);
    }
    else if (fmpz_mpoly_is_one(t2, ctx))
    {
        fmpz_zero(mu1);
        fmpz_one(mu2);

        fmpz_mpoly_swap(hbar, bbar, ctx);

        fmpz_mpolyl_lead_coeff(Hlcs + n, B, 1, ctx);
        fmpz_mpolyl_lead_coeff(t2, A, 1, ctx);
        success = fmpz_mpoly_gcd(Glcs + n, Hlcs + n, t2, ctx);
        if (!success)
            goto cleanup;

        H = (fmpz_mpoly_struct *) B;

        gamma_is_one = fmpz_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                fmpz_mpoly_swap(Hevals + i, Bevals + i, ctx);
    }
    else
    {
        int mu_tries_remaining = 10;

    next_mu:

        if (--mu_tries_remaining < 0)
        {
            success = 0;
            goto cleanup;
        }

        fmpq_next_signed_calkin_wilf(mu, mu);
        fmpz_set(mu1, fmpq_numref(mu));
        fmpz_set(mu2, fmpq_denref(mu));

        fmpz_mpoly_scalar_fmma(hbar, abar, mu1, bbar, mu2, ctx);

        /* make sure the linear combo did not drop degree */
        if (fmpz_mpoly_degree_si(hbar, 0, ctx) != FLINT_MAX(Adegx, Bdegx) - gdegx)
            goto next_mu;

        /* make sure the linear combo is prime to g */
        success = fmpz_mpoly_gcd(t1, hbar, g, ctx);
        if (!success)
            goto cleanup;

        if (!fmpz_mpoly_is_fmpz(t1, ctx))
            goto next_mu;

        fmpz_mpolyl_lead_coeff(t1, A, 1, ctx);
        fmpz_mpolyl_lead_coeff(t2, B, 1, ctx);
        success = fmpz_mpoly_gcd(Glcs + n, t1, t2, ctx);
        if (!success)
            goto cleanup;

        H = Hevals + n;
        fmpz_mpoly_scalar_fmma(H, A, mu1, B, mu2, ctx);
        fmpz_mpolyl_lead_coeff(Hlcs + n, H, 1, ctx);

        gamma_is_one = fmpz_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                fmpz_mpoly_scalar_fmma(Hevals + i, Aevals + i, mu1,
                                                   Bevals + i, mu2, ctx);
    }

    if (!gamma_is_one)
    {
        fmpz_mpoly_mul(Hevals + n, H, Glcs + n, ctx);
        H = Hevals + n;
        for (i = n - 1; i >= 0; i--)
            fmpz_mpoly_evaluate_one_fmpz(Hevals + i, Hevals + i + 1,
                                                       i + 1, alphas + i, ctx);
    }

    success = H->bits <= FLINT_BITS ||
              fmpz_mpoly_repack_bits_inplace(H, FLINT_BITS, ctx);
    if (!success)
        goto cleanup;

    /* the evals should all fit in H->bits */
    for (i = 0; i < n; i++)
        fmpz_mpoly_repack_bits_inplace(Hevals + i, H->bits, ctx);

    fmpz_mpoly_degrees_si(Hdegs, H, ctx);

    /* computed evaluated leading coeffs */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Glcs + i, Glcs + i + 1, i + 1, alphas + i, ctx);
        fmpz_mpoly_evaluate_one_fmpz(Hlcs + i, Hlcs + i + 1, i + 1, alphas + i, ctx);
        /* evaluation could have killed gamma */
        if (fmpz_mpoly_is_zero(Glcs + i, ctx) ||
            fmpz_mpoly_is_zero(Hlcs + i, ctx))
        {
            goto next_alpha;
        }
    }

    /* make the lcs match Glcs[0], Hlcs[0], g & hbar not nec primitive */

    fmpz_mpoly_scalar_mul_fmpz(Hfac + 0, g, Glcs[0].coeffs + 0, ctx);
    success = fmpz_mpoly_scalar_divides_fmpz(Hfac + 0, Hfac + 0, g->coeffs + 0, ctx);
    FLINT_ASSERT(success);

    fmpz_mpoly_scalar_mul_fmpz(Hfac + 1, hbar, Hlcs[0].coeffs + 0, ctx);
    success = fmpz_mpoly_scalar_divides_fmpz(Hfac + 1, Hfac + 1, hbar->coeffs + 0, ctx);
    FLINT_ASSERT(success);

    for (k = 1; k <= n; k++)
    {
        _fmpz_mpoly_set_lead0(Htfac + 0, Hfac + 0, Glcs + k, ctx);
        _fmpz_mpoly_set_lead0(Htfac + 1, Hfac + 1, Hlcs + k, ctx);
        success = fmpz_mpoly_hlift(k, Htfac, 2, alphas,
                                           k < n ? Hevals + k : H, Hdegs, ctx);
        if (!success)
            goto next_alpha;

        fmpz_mpoly_swap(Hfac + 0, Htfac + 0, ctx);
        fmpz_mpoly_swap(Hfac + 1, Htfac + 1, ctx);
    }

    success = fmpz_mpolyl_content(t1, Hfac + 0, 1, ctx);
    if (!success)
        goto cleanup;

    success = fmpz_mpoly_divides(G, Hfac + 0, t1, ctx);
    FLINT_ASSERT(success);

    if (fmpz_is_zero(mu2))
    {
        FLINT_ASSERT(fmpz_is_one(mu1));
        /* the division by t1 should succeed, but let's be careful */
        fmpz_mpolyl_lead_coeff(t1, G, 1, ctx);
        success = fmpz_mpoly_divides(Abar, Hfac + 1, t1, ctx) &&
                  fmpz_mpoly_divides(Bbar, B, G, ctx);
    }
    else if (fmpz_is_zero(mu1))
    {
        FLINT_ASSERT(fmpz_is_one(mu2));
        /* ditto */
        fmpz_mpolyl_lead_coeff(t1, G, 1, ctx);
        success = fmpz_mpoly_divides(Bbar, Hfac + 1, t1, ctx) &&
                  fmpz_mpoly_divides(Abar, A, G, ctx);
    }
    else
    {
        success = fmpz_mpoly_divides(Abar, A, G, ctx) &&
                  fmpz_mpoly_divides(Bbar, B, G, ctx);
    }

    if (!success)
        goto next_alpha;

    success = 1;

cleanup:

    flint_randclear(state);

    flint_free(Hdegs);

    for (i = 0; i < n + 1; i++)
    {
        fmpz_mpoly_clear(Glcs + i, ctx);
        fmpz_mpoly_clear(Hlcs + i, ctx);
        fmpz_mpoly_clear(Hevals + i, ctx);
    }
    flint_free(Glcs);

    _fmpz_vec_clear(alphas, 2*n);

    for (i = 0; i < n; i++)
    {
        fmpz_mpoly_clear(Aevals + i, ctx);
        fmpz_mpoly_clear(Bevals + i, ctx);
    }
    flint_free(Aevals);

    fmpz_clear(q);
    fmpq_clear(mu);
    fmpz_clear(mu1);
    fmpz_clear(mu2);

    fmpz_mpoly_clear(t1, ctx);
    fmpz_mpoly_clear(t2, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(abar, ctx);
    fmpz_mpoly_clear(bbar, ctx);
    fmpz_mpoly_clear(hbar, ctx);

    fmpz_mpoly_clear(Hfac + 0, ctx);
    fmpz_mpoly_clear(Hfac + 1, ctx);
    fmpz_mpoly_clear(Htfac + 0, ctx);
    fmpz_mpoly_clear(Htfac + 1, ctx);

    if (success)
    {
        fmpz_mpoly_repack_bits_inplace(G, bits, ctx);
        fmpz_mpoly_repack_bits_inplace(Abar, bits, ctx);
        fmpz_mpoly_repack_bits_inplace(Bbar, bits, ctx);

        FLINT_ASSERT(G->length > 0);
        FLINT_ASSERT(Abar->length > 0);
        FLINT_ASSERT(Bbar->length > 0);
    }

    return success;
}

