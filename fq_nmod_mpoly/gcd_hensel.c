/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"
#include "fq_nmod_mpoly_factor.h"


/*
    in F[x_1, ..., x_n][X]:

    cont_X(A) = cont_X(B) = 1
    compute gamma = gcd(lc_X(A), lc_X(B))

    try to find evaluation point x_i -> alpha_i
    a = A(x_i = alpha_i) in F[X]
    b = B(x_i = alpha_i) in F[X]

    compute univariate a = g*abar, b = g*bbar

    try to find mu1, mu2 in F with gcd(g, mu1*abar + mu2*bbar) = 1

    set H = mu1*A + mu2*B
    lift the univariate factorization g * (mu1*abar + mu2*bbar) against gamma*H
    imposing (gamma, lc_X H) as the leading coefficients

    remove content from the lift of g to get G and divide to get Abar and Bbar
*/
int fq_nmod_mpolyl_gcd_hensel_smprime(
    fq_nmod_mpoly_t G, slong Gdeg,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success, alphas_tries_remaining, gamma_is_one;
    const slong n = ctx->minfo->nvars - 1;
    slong i, k;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_struct * alphas;
    fq_nmod_t mu1, mu2;
    fq_nmod_mpoly_struct * Aevals, * Bevals, * Hevals;
    fq_nmod_mpoly_struct * H; /* points to A, B, or Hevals + n */
    fq_nmod_mpoly_struct * Glcs, * Hlcs;
    fq_nmod_mpoly_struct Hfac[2], Htfac[2];
    slong * Hdegs;
    slong Adegx, Bdegx, gdegx;
    fq_nmod_mpoly_t t1, t2, g, abar, bbar, hbar;
    flint_rand_t state;
    mp_limb_t * tmp, * q;

    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    flint_randinit(state);

    tmp = FLINT_ARRAY_ALLOC(d*(N_FQ_MUL_INV_ITCH + 1), mp_limb_t);
    q = tmp + d*N_FQ_MUL_INV_ITCH;

    fq_nmod_init(mu1, ctx->fqctx);
    fq_nmod_init(mu2, ctx->fqctx);

    Hdegs  = FLINT_ARRAY_ALLOC(n + 1, slong);

    Glcs   = FLINT_ARRAY_ALLOC(3*(n + 1), fq_nmod_mpoly_struct);
    Hlcs   = Glcs + (n + 1);
    Hevals = Hlcs + (n + 1);
    for (i = 0; i < n + 1; i++)
    {
        fq_nmod_mpoly_init(Glcs + i, ctx);
        fq_nmod_mpoly_init(Hlcs + i, ctx);
        fq_nmod_mpoly_init(Hevals + i, ctx);
    }

	alphas = FLINT_ARRAY_ALLOC(n, fq_nmod_struct);
    Aevals = FLINT_ARRAY_ALLOC(2*(n + 1), fq_nmod_mpoly_struct);
    Bevals = Aevals + (n + 1);
	for (i = 0; i < n; i++)
    {
        fq_nmod_init(alphas + i, ctx->fqctx);
		fq_nmod_mpoly_init(Aevals + i, ctx);
		fq_nmod_mpoly_init(Bevals + i, ctx);
    }

    fq_nmod_mpoly_init(t1, ctx);
    fq_nmod_mpoly_init(t2, ctx);
    fq_nmod_mpoly_init(g, ctx);
    fq_nmod_mpoly_init(abar, ctx);
    fq_nmod_mpoly_init(bbar, ctx);
    fq_nmod_mpoly_init(hbar, ctx);

    fq_nmod_mpoly_init(Hfac + 0, ctx);
    fq_nmod_mpoly_init(Hfac + 1, ctx);
    fq_nmod_mpoly_init(Htfac + 0, ctx);
    fq_nmod_mpoly_init(Htfac + 1, ctx);

    /* init done */

    alphas_tries_remaining = 10;

    /* try all zeros first */
    for (i = 0; i < n; i++)
        fq_nmod_zero(alphas + i, ctx->fqctx);

    goto got_alpha;

next_alpha:

    if (--alphas_tries_remaining < 0)
	{
		success = 0;
        goto cleanup;
	}

    for (i = 0; i < n; i++)
        fq_nmod_rand(alphas + i, state, ctx->fqctx);

got_alpha:

    /* ensure deg_X do not drop under evaluation */
    Adegx = fq_nmod_mpoly_degree_si(A, 0, ctx);
    Bdegx = fq_nmod_mpoly_degree_si(B, 0, ctx);

	for (i = n - 1; i >= 0; i--)
	{
		fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i, i == n - 1 ? A :
                                       Aevals + i + 1, i + 1, alphas + i, ctx);
		fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + i, i == n - 1 ? B :
                                       Bevals + i + 1, i + 1, alphas + i, ctx);
		if (Adegx != fq_nmod_mpoly_degree_si(Aevals + i, 0, ctx) ||
            Bdegx != fq_nmod_mpoly_degree_si(Bevals + i, 0, ctx))
        {
    		goto next_alpha;
        }
	}

    /* univariate gcd */
    success = fq_nmod_mpoly_gcd_cofactors(g, abar, bbar,
                                                Aevals + 0, Bevals + 0, ctx) &&
              fq_nmod_mpoly_gcd(t1, g, abar, ctx) &&
              fq_nmod_mpoly_gcd(t2, g, bbar, ctx);
    if (!success)
        goto cleanup;

    gdegx = fq_nmod_mpoly_degree_si(g, 0, ctx);

    if (gdegx == 0)
    {
        fq_nmod_mpoly_set(Abar, A, ctx);
        fq_nmod_mpoly_set(Bbar, B, ctx);
        fq_nmod_mpoly_one(G, ctx);
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
        goto next_alpha;
    }

    /* the degbound gdegx (== Gdeg) has at least two witnesses now */

    if (gdegx == Adegx)
    {
        if (fq_nmod_mpoly_divides(Bbar, B, A, ctx))
        {
            fq_nmod_mpoly_set(G, A, ctx);
            fq_nmod_mpoly_one(Abar, ctx);
            success = 1;
            goto cleanup;
        }

        goto next_alpha;
    }
    else if (gdegx == Bdegx)
    {
        if (fq_nmod_mpoly_divides(Abar, A, B, ctx))
        {
            fq_nmod_mpoly_set(G, B, ctx);
            fq_nmod_mpoly_one(Bbar, ctx);
            success = 1;
            goto cleanup;
        }

        goto next_alpha;
    }

    FLINT_ASSERT(0 < gdegx && gdegx < FLINT_MIN(Adegx, Bdegx));

    /* set Hlcs[n], Glcs[n] (gamma), H, and Hevals */
    if (fq_nmod_mpoly_is_one(t1, ctx))
    {
        fq_nmod_one(mu1, ctx->fqctx);
        fq_nmod_zero(mu2, ctx->fqctx);

        fq_nmod_mpoly_swap(hbar, abar, ctx);

        fq_nmod_mpolyl_lead_coeff(Hlcs + n, A, 1, ctx);
        fq_nmod_mpolyl_lead_coeff(t2, B, 1, ctx);
        success = fq_nmod_mpoly_gcd(Glcs + n, Hlcs + n, t2, ctx);
        if (!success)
            goto cleanup;

        H = (fq_nmod_mpoly_struct *) A;

        gamma_is_one = fq_nmod_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                fq_nmod_mpoly_swap(Hevals + i, Aevals + i, ctx);
    }
    else if (fq_nmod_mpoly_is_one(t2, ctx))
    {
        fq_nmod_zero(mu1, ctx->fqctx);
        fq_nmod_one(mu2, ctx->fqctx);

        fq_nmod_mpoly_swap(hbar, bbar, ctx);

        fq_nmod_mpolyl_lead_coeff(Hlcs + n, B, 1, ctx);
        fq_nmod_mpolyl_lead_coeff(t2, A, 1, ctx);
        success = fq_nmod_mpoly_gcd(Glcs + n, Hlcs + n, t2, ctx);
        if (!success)
            goto cleanup;

        H = (fq_nmod_mpoly_struct *) B;

        gamma_is_one = fq_nmod_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                fq_nmod_mpoly_swap(Hevals + i, Bevals + i, ctx);
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

        fq_nmod_one(mu1, ctx->fqctx);
        fq_nmod_rand_not_zero(mu2, state, ctx->fqctx);
        fq_nmod_mpoly_scalar_addmul_fq_nmod(hbar, abar, bbar, mu2, ctx);

        /* make sure the linear combo did not drop degree */
        if (fq_nmod_mpoly_degree_si(hbar, 0, ctx) != FLINT_MAX(Adegx, Bdegx) - gdegx)
            goto next_mu;

        /* make sure the linear combo is prime to g */
        success = fq_nmod_mpoly_gcd(t1, hbar, g, ctx);
        if (!success)
            goto cleanup;

        if (!fq_nmod_mpoly_is_one(t1, ctx))
            goto next_mu;

        fq_nmod_mpolyl_lead_coeff(t1, A, 1, ctx);
        fq_nmod_mpolyl_lead_coeff(t2, B, 1, ctx);
        success = fq_nmod_mpoly_gcd(Glcs + n, t1, t2, ctx);
        if (!success)
            goto cleanup;

        H = Hevals + n;
        fq_nmod_mpoly_scalar_addmul_fq_nmod(H, A, B, mu2, ctx);
        fq_nmod_mpolyl_lead_coeff(Hlcs + n, H, 1, ctx);

        gamma_is_one = fq_nmod_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                fq_nmod_mpoly_scalar_addmul_fq_nmod(Hevals + i, Aevals + i,
                                                         Bevals + i, mu2, ctx);
    }

    if (!gamma_is_one)
    {
        fq_nmod_mpoly_mul(Hevals + n, H, Glcs + n, ctx);
        H = Hevals + n;
        for (i = n - 1; i >= 0; i--)
            fq_nmod_mpoly_evaluate_one_fq_nmod(Hevals + i, Hevals + i + 1,
                                                       i + 1, alphas + i, ctx);
    }

    success = H->bits <= FLINT_BITS ||
              fq_nmod_mpoly_repack_bits_inplace(H, FLINT_BITS, ctx);
    if (!success)
        goto cleanup;

    /* the evals should all fit in H->bits */
    for (i = 0; i < n; i++)
        fq_nmod_mpoly_repack_bits_inplace(Hevals + i, H->bits, ctx);

    fq_nmod_mpoly_degrees_si(Hdegs, H, ctx);

    /* computed evaluated leading coeffs */
    for (i = n - 1; i >= 0; i--)
    {
        fq_nmod_mpoly_evaluate_one_fq_nmod(Glcs + i, Glcs + i + 1,
                                                       i + 1, alphas + i, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Hlcs + i, Hlcs + i + 1,
                                                       i + 1, alphas + i, ctx);
        /* evaluation could have killed gamma */
        if (fq_nmod_mpoly_is_zero(Glcs + i, ctx) ||
            fq_nmod_mpoly_is_zero(Hlcs + i, ctx))
        {
            goto next_alpha;
        }
    }

    /* make the leading coefficients match Glcs[0], Hlcs[0] */
    FLINT_ASSERT(fq_nmod_mpoly_is_fq_nmod(Glcs + 0, ctx) && Glcs[0].length == 1);
    FLINT_ASSERT(fq_nmod_mpoly_is_fq_nmod(Hlcs + 0, ctx) && Hlcs[0].length == 1);

    _n_fq_inv(q, g->coeffs + 0, ctx->fqctx, tmp);
    _n_fq_mul(q, q, Glcs[0].coeffs + 0, ctx->fqctx, tmp);
    fq_nmod_mpoly_scalar_mul_n_fq(Hfac + 0, g, q, ctx);

    _n_fq_inv(q, hbar->coeffs + 0, ctx->fqctx, tmp);
    _n_fq_mul(q, q, Hlcs[0].coeffs + 0, ctx->fqctx, tmp);
    fq_nmod_mpoly_scalar_mul_n_fq(Hfac + 1, hbar, q, ctx);

    for (k = 1; k <= n; k++)
    {
        _fq_nmod_mpoly_set_lead0(Htfac + 0, Hfac + 0, Glcs + k, ctx);
        _fq_nmod_mpoly_set_lead0(Htfac + 1, Hfac + 1, Hlcs + k, ctx);
        success = fq_nmod_mpoly_hlift(k, Htfac, 2, alphas,
                                           k < n ? Hevals + k : H, Hdegs, ctx);
        if (!success)
            goto next_alpha;

        fq_nmod_mpoly_swap(Hfac + 0, Htfac + 0, ctx);
        fq_nmod_mpoly_swap(Hfac + 1, Htfac + 1, ctx);
    }

    success = fq_nmod_mpolyl_content(t1, Hfac + 0, 1, ctx);
    if (!success)
        goto cleanup;

    success = fq_nmod_mpoly_divides(G, Hfac + 0, t1, ctx);
    FLINT_ASSERT(success);

    if (fq_nmod_is_zero(mu2, ctx->fqctx))
    {
        FLINT_ASSERT(fq_nmod_is_one(mu1, ctx->fqctx));
        /* the division by t1 should succeed, but let's be careful */
        fq_nmod_mpolyl_lead_coeff(t1, G, 1, ctx);
        success = fq_nmod_mpoly_divides(Abar, Hfac + 1, t1, ctx) &&
                  fq_nmod_mpoly_divides(Bbar, B, G, ctx);
    }
    else if (fq_nmod_is_zero(mu1, ctx->fqctx))
    {
        FLINT_ASSERT(fq_nmod_is_one(mu2, ctx->fqctx));
        /* ditto */
        fq_nmod_mpolyl_lead_coeff(t1, G, 1, ctx);
        success = fq_nmod_mpoly_divides(Bbar, Hfac + 1, t1, ctx) &&
                  fq_nmod_mpoly_divides(Abar, A, G, ctx);
    }
    else
    {
        FLINT_ASSERT(fq_nmod_is_one(mu1, ctx->fqctx));
        success = fq_nmod_mpoly_divides(Abar, A, G, ctx) &&
                  fq_nmod_mpoly_divides(Bbar, B, G, ctx);
    }

    if (!success)
        goto next_alpha;

    success = 1;

cleanup:

    flint_randclear(state);

    flint_free(tmp);

    fq_nmod_clear(mu1, ctx->fqctx);
    fq_nmod_clear(mu2, ctx->fqctx);

    flint_free(Hdegs);

    for (i = 0; i < n + 1; i++)
    {
        fq_nmod_mpoly_clear(Glcs + i, ctx);
        fq_nmod_mpoly_clear(Hlcs + i, ctx);
        fq_nmod_mpoly_clear(Hevals + i, ctx);
    }
    flint_free(Glcs);

	for (i = 0; i < n; i++)
    {
        fq_nmod_clear(alphas + i, ctx->fqctx);
		fq_nmod_mpoly_clear(Aevals + i, ctx);
		fq_nmod_mpoly_clear(Bevals + i, ctx);
    }
    flint_free(alphas);
    flint_free(Aevals);

    fq_nmod_mpoly_clear(t1, ctx);
    fq_nmod_mpoly_clear(t2, ctx);
    fq_nmod_mpoly_clear(g, ctx);
    fq_nmod_mpoly_clear(abar, ctx);
    fq_nmod_mpoly_clear(bbar, ctx);
    fq_nmod_mpoly_clear(hbar, ctx);

    fq_nmod_mpoly_clear(Hfac + 0, ctx);
    fq_nmod_mpoly_clear(Hfac + 1, ctx);
    fq_nmod_mpoly_clear(Htfac + 0, ctx);
    fq_nmod_mpoly_clear(Htfac + 1, ctx);

    if (success)
    {
        fq_nmod_mpoly_repack_bits_inplace(G, bits, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
    }

	return success;
}


/* should find its way back here in interesting cases */
int fq_nmod_mpoly_gcd_hensel(
    fq_nmod_mpoly_t G,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    if (fq_nmod_mpoly_is_zero(A, ctx) || fq_nmod_mpoly_is_zero(B, ctx))
        return fq_nmod_mpoly_gcd(G, A, B, ctx);

    return _fq_nmod_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_HENSEL);
}

