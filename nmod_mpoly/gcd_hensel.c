/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "fq_zech_mpoly.h"
#include "fq_zech_mpoly_factor.h"

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
int nmod_mpolyl_gcd_hensel_smprime(
    nmod_mpoly_t G, slong Gdeg, /* upperbound on deg_X(G) */
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int success, alphas_tries_remaining, gamma_is_one;
    const slong n = ctx->minfo->nvars - 1;
    slong i, k;
    flint_bitcnt_t bits = A->bits;
    mp_limb_t * alphas, * prev_alphas;
    mp_limb_t q, mu1, mu2;
    nmod_mpoly_struct * Aevals, * Bevals, * Hevals;
    nmod_mpoly_struct * H; /* points to A, B, or Hevals + n */
    nmod_mpoly_struct * Glcs, * Hlcs;
    nmod_mpoly_struct Hfac[2], Htfac[2];
    slong * Hdegs;
    slong Adegx, Bdegx, gdegx;
    nmod_mpoly_t t1, t2, g, abar, bbar, hbar;
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

    Glcs   = FLINT_ARRAY_ALLOC(3*(n + 1), nmod_mpoly_struct);
    Hlcs   = Glcs + (n + 1);
    Hevals = Hlcs + (n + 1);
    for (i = 0; i < n + 1; i++)
    {
        nmod_mpoly_init(Glcs + i, ctx);
        nmod_mpoly_init(Hlcs + i, ctx);
        nmod_mpoly_init(Hevals + i, ctx);
    }

	alphas = FLINT_ARRAY_ALLOC(2*n, mp_limb_t);
    prev_alphas = alphas + n;
    Aevals = FLINT_ARRAY_ALLOC(2*(n + 1), nmod_mpoly_struct);
    Bevals = Aevals + (n + 1);
	for (i = 0; i < n; i++)
    {
		nmod_mpoly_init(Aevals + i, ctx);
		nmod_mpoly_init(Bevals + i, ctx);
    }

    nmod_mpoly_init(t1, ctx);
    nmod_mpoly_init(t2, ctx);
    nmod_mpoly_init(g, ctx);
    nmod_mpoly_init(abar, ctx);
    nmod_mpoly_init(bbar, ctx);
    nmod_mpoly_init(hbar, ctx);

    nmod_mpoly_init(Hfac + 0, ctx);
    nmod_mpoly_init(Hfac + 1, ctx);
    nmod_mpoly_init(Htfac + 0, ctx);
    nmod_mpoly_init(Htfac + 1, ctx);

    /* init done */

    alphas_tries_remaining = 10;

    /* try all zeros first */
    for (i = 0; i < n; i++)
    {
        prev_alphas[i] = ctx->mod.n; /* no previous at this point */
        alphas[i] = 0;
    }

    goto got_alpha;

next_alpha:

    if (--alphas_tries_remaining < 0)
	{
		success = 0;
        goto cleanup;
	}

    for (i = 0; i < n; i++)
    {
        do {
            alphas[i] = n_urandint(state, ctx->mod.n);
        } while (alphas[i] == prev_alphas[i]);
    }

got_alpha:

    /* ensure deg_X do not drop under evaluation */
    Adegx = nmod_mpoly_degree_si(A, 0, ctx);
    Bdegx = nmod_mpoly_degree_si(B, 0, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        nmod_mpoly_evaluate_one_ui(Aevals + i, i == n - 1 ? A :
                                        Aevals + i + 1, i + 1, alphas[i], ctx);
        nmod_mpoly_evaluate_one_ui(Bevals + i, i == n - 1 ? B :
                                        Bevals + i + 1, i + 1, alphas[i], ctx);
        if (Adegx != nmod_mpoly_degree_si(Aevals + i, 0, ctx) ||
            Bdegx != nmod_mpoly_degree_si(Bevals + i, 0, ctx))
        {
            goto next_alpha;
        }
	}

    /* univariate gcd */
    success = nmod_mpoly_gcd_cofactors(g, abar, bbar,
                                       Aevals + 0, Bevals + 0, ctx) &&
              nmod_mpoly_gcd(t1, g, abar, ctx) &&
              nmod_mpoly_gcd(t2, g, bbar, ctx);
    if (!success)
        goto cleanup;

    gdegx = nmod_mpoly_degree_si(g, 0, ctx);

    if (gdegx == 0)
    {
        /* G is trivial */
        nmod_mpoly_set(Abar, A, ctx);
        nmod_mpoly_set(Bbar, B, ctx);
        nmod_mpoly_one(G, ctx);
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
            prev_alphas[i] = alphas[i];

        goto next_alpha;
    }

    /* the degbound gdegx (== Gdeg) has at least two witnesses now */

    if (gdegx == Adegx)
    {
        if (nmod_mpoly_divides(Bbar, B, A, ctx))
        {
            nmod_mpoly_set(G, A, ctx);
            nmod_mpoly_one(Abar, ctx);
            success = 1;
            goto cleanup;
        }

        goto next_alpha;
    }
    else if (gdegx == Bdegx)
    {
        if (nmod_mpoly_divides(Abar, A, B, ctx))
        {
            nmod_mpoly_set(G, B, ctx);
            nmod_mpoly_one(Bbar, ctx);
            success = 1;
            goto cleanup;
        }

        goto next_alpha;
    }

    FLINT_ASSERT(0 < gdegx && gdegx < FLINT_MIN(Adegx, Bdegx));

    /* set Hlcs[n], Glcs[n] (gamma), H, and Hevals */
    if (nmod_mpoly_is_one(t1, ctx))
    {
        mu1 = 1;
        mu2 = 0;

        nmod_mpoly_swap(hbar, abar, ctx);

        nmod_mpolyl_lead_coeff(Hlcs + n, A, 1, ctx);
        nmod_mpolyl_lead_coeff(t2, B, 1, ctx);
        success = nmod_mpoly_gcd(Glcs + n, Hlcs + n, t2, ctx);
        if (!success)
            goto cleanup;

        H = (nmod_mpoly_struct *) A;

        gamma_is_one = nmod_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                nmod_mpoly_swap(Hevals + i, Aevals + i, ctx);
    }
    else if (nmod_mpoly_is_one(t2, ctx))
    {
        mu1 = 0;
        mu2 = 1;

        nmod_mpoly_swap(hbar, bbar, ctx);

        nmod_mpolyl_lead_coeff(Hlcs + n, B, 1, ctx);
        nmod_mpolyl_lead_coeff(t2, A, 1, ctx);
        success = nmod_mpoly_gcd(Glcs + n, Hlcs + n, t2, ctx);
        if (!success)
            goto cleanup;

        H = (nmod_mpoly_struct *) B;

        gamma_is_one = nmod_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                nmod_mpoly_swap(Hevals + i, Bevals + i, ctx);
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

        mu1 = 1;
        mu2 = n_urandint(state, ctx->mod.n - 1) + 1;
        nmod_mpoly_scalar_addmul_ui(hbar, abar, bbar, mu2, ctx);

        /* make sure the linear combo did not drop degree */
        if (nmod_mpoly_degree_si(hbar, 0, ctx) != FLINT_MAX(Adegx, Bdegx) - gdegx)
            goto next_mu;

        /* make sure the linear combo is prime to g */
        success = nmod_mpoly_gcd(t1, hbar, g, ctx);
        if (!success)
            goto cleanup;

        if (!nmod_mpoly_is_one(t1, ctx))
            goto next_mu;

        nmod_mpolyl_lead_coeff(t1, A, 1, ctx);
        nmod_mpolyl_lead_coeff(t2, B, 1, ctx);
        success = nmod_mpoly_gcd(Glcs + n, t1, t2, ctx);
        if (!success)
            goto cleanup;

        H = Hevals + n;
        nmod_mpoly_scalar_addmul_ui(H, A, B, mu2, ctx);
        nmod_mpolyl_lead_coeff(Hlcs + n, H, 1, ctx);

        gamma_is_one = nmod_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                nmod_mpoly_scalar_addmul_ui(Hevals + i, Aevals + i,
                                                        Bevals + i, mu2, ctx);
    }

    if (!gamma_is_one)
    {
        nmod_mpoly_mul(Hevals + n, H, Glcs + n, ctx);
        H = Hevals + n;
        for (i = n - 1; i >= 0; i--)
            nmod_mpoly_evaluate_one_ui(Hevals + i, Hevals + i + 1,
                                                        i + 1, alphas[i], ctx);
    }

    success = H->bits <= FLINT_BITS ||
              nmod_mpoly_repack_bits_inplace(H, FLINT_BITS, ctx);
    if (!success)
        goto cleanup;

    /* the evals should all fit in H->bits */
    for (i = 0; i < n; i++)
        nmod_mpoly_repack_bits_inplace(Hevals + i, H->bits, ctx);

    nmod_mpoly_degrees_si(Hdegs, H, ctx);

    /* computed evaluated leading coeffs */
    for (i = n - 1; i >= 0; i--)
    {
        nmod_mpoly_evaluate_one_ui(Glcs + i, Glcs + i + 1, i + 1, alphas[i], ctx);
        nmod_mpoly_evaluate_one_ui(Hlcs + i, Hlcs + i + 1, i + 1, alphas[i], ctx);
        /* evaluation could have killed gamma */
        if (nmod_mpoly_is_zero(Glcs + i, ctx) ||
            nmod_mpoly_is_zero(Hlcs + i, ctx))
        {
            goto next_alpha;
        }
    }

    /* make the leading coefficients match Glcs[0], Hlcs[0] */

    FLINT_ASSERT(nmod_mpoly_is_ui(Glcs + 0, ctx) && Glcs[0].length == 1);
    FLINT_ASSERT(nmod_mpoly_is_ui(Hlcs + 0, ctx) && Hlcs[0].length == 1);

    q = nmod_inv(g->coeffs[0], ctx->mod);
    q = nmod_mul(q, Glcs[0].coeffs[0], ctx->mod);
    nmod_mpoly_scalar_mul_nmod_invertible(Hfac + 0, g, q, ctx);

    q = nmod_inv(hbar->coeffs[0], ctx->mod);
    q = nmod_mul(q, Hlcs[0].coeffs[0], ctx->mod);
    nmod_mpoly_scalar_mul_nmod_invertible(Hfac + 1, hbar, q, ctx);

    for (k = 1; k <= n; k++)
    {
        _nmod_mpoly_set_lead0(Htfac + 0, Hfac + 0, Glcs + k, ctx);
        _nmod_mpoly_set_lead0(Htfac + 1, Hfac + 1, Hlcs + k, ctx);
        success = nmod_mpoly_hlift(k, Htfac, 2, alphas,
                                           k < n ? Hevals + k : H, Hdegs, ctx);
        if (!success)
            goto next_alpha;

        nmod_mpoly_swap(Hfac + 0, Htfac + 0, ctx);
        nmod_mpoly_swap(Hfac + 1, Htfac + 1, ctx);
    }

    success = nmod_mpolyl_content(t1, Hfac + 0, 1, ctx);
    if (!success)
        goto cleanup;

    success = nmod_mpoly_divides(G, Hfac + 0, t1, ctx);
    FLINT_ASSERT(success);

    if (mu2 == 0)
    {
        FLINT_ASSERT(mu1 == 1);
        /* the division by t1 should succeed, but let's be careful */
        nmod_mpolyl_lead_coeff(t1, G, 1, ctx);
        success = nmod_mpoly_divides(Abar, Hfac + 1, t1, ctx) &&
                  nmod_mpoly_divides(Bbar, B, G, ctx);
    }
    else if (mu1 == 0)
    {
        FLINT_ASSERT(mu2 == 1);
        /* ditto */
        nmod_mpolyl_lead_coeff(t1, G, 1, ctx);
        success = nmod_mpoly_divides(Bbar, Hfac + 1, t1, ctx) &&
                  nmod_mpoly_divides(Abar, A, G, ctx);
    }
    else
    {
        FLINT_ASSERT(mu1 == 1);
        success = nmod_mpoly_divides(Abar, A, G, ctx) &&
                  nmod_mpoly_divides(Bbar, B, G, ctx);
    }

    if (!success)
        goto next_alpha;

    success = 1;

cleanup:

    flint_randclear(state);

    flint_free(Hdegs);

    for (i = 0; i < n + 1; i++)
    {
        nmod_mpoly_clear(Glcs + i, ctx);
        nmod_mpoly_clear(Hlcs + i, ctx);
        nmod_mpoly_clear(Hevals + i, ctx);
    }
    flint_free(Glcs);

	for (i = 0; i < n; i++)
    {
		nmod_mpoly_clear(Aevals + i, ctx);
		nmod_mpoly_clear(Bevals + i, ctx);
    }
    flint_free(alphas);
    flint_free(Aevals);

    nmod_mpoly_clear(t1, ctx);
    nmod_mpoly_clear(t2, ctx);
    nmod_mpoly_clear(g, ctx);
    nmod_mpoly_clear(abar, ctx);
    nmod_mpoly_clear(bbar, ctx);
    nmod_mpoly_clear(hbar, ctx);

    nmod_mpoly_clear(Hfac + 0, ctx);
    nmod_mpoly_clear(Hfac + 1, ctx);
    nmod_mpoly_clear(Htfac + 0, ctx);
    nmod_mpoly_clear(Htfac + 1, ctx);

    if (success)
    {
        nmod_mpoly_repack_bits_inplace(G, bits, ctx);
        nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
        nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);

        FLINT_ASSERT(G->length > 0);
        FLINT_ASSERT(Abar->length > 0);
        FLINT_ASSERT(Bbar->length > 0);
    }

	return success;
}


int nmod_mpolyl_gcd_hensel_medprime(
    nmod_mpoly_t G, slong Gdeg, /* upperbound on deg_X(G) */
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t smA,
    const nmod_mpoly_t smB,
    const nmod_mpoly_ctx_t smctx)
{
    int success, alphas_tries_remaining, gamma_is_one;
    const slong n = smctx->minfo->nvars - 1;
    slong i, k;
    flint_bitcnt_t bits = smA->bits;
    fq_zech_struct * alphas;
    fq_zech_t q, mu1, mu2;
    fq_zech_mpoly_t A, B;
    fq_zech_mpoly_struct * Aevals, * Bevals, * Hevals;
    fq_zech_mpoly_struct * H; /* points to A, B, or Hevals + n */
    fq_zech_mpoly_struct * Glcs, * Hlcs;
    fq_zech_mpoly_struct Hfac[2], Htfac[2];
    nmod_mpoly_struct smHfac[2];
    slong * Hdegs;
    slong Adegx, Bdegx, gdegx;
    fq_zech_mpoly_t t1, t2, g, abar, bbar, hbar;
    nmod_mpoly_t t;
    flint_rand_t state;
    fq_zech_mpoly_ctx_t ctx;
    slong edeg, max_degree = n_flog(1000000, smctx->mod.n);

    if (max_degree < 2)
        return 0;

    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(smA->length > 0);
    FLINT_ASSERT(smB->length > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(smA->bits == bits);
    FLINT_ASSERT(smB->bits == bits);
    FLINT_ASSERT(smctx->minfo->ord == ORD_LEX);

    edeg = 1 + n_clog(500, smctx->mod.n);
    edeg = FLINT_MAX(2, edeg);
    edeg = FLINT_MIN(edeg, max_degree);
    fq_zech_mpoly_ctx_init_deg(ctx, smctx->minfo->nvars, ORD_LEX, smctx->mod.n, edeg);

    fq_zech_init(q, ctx->fqctx);
    fq_zech_init(mu1, ctx->fqctx);
    fq_zech_init(mu2, ctx->fqctx);
    fq_zech_mpoly_init(A, ctx);
    fq_zech_mpoly_init(B, ctx);

    flint_randinit(state);

    Hdegs  = FLINT_ARRAY_ALLOC(n + 1, slong);

    Glcs   = FLINT_ARRAY_ALLOC(3*(n + 1), fq_zech_mpoly_struct);
    Hlcs   = Glcs + (n + 1);
    Hevals = Hlcs + (n + 1);
    for (i = 0; i < n + 1; i++)
    {
        fq_zech_mpoly_init(Glcs + i, ctx);
        fq_zech_mpoly_init(Hlcs + i, ctx);
        fq_zech_mpoly_init(Hevals + i, ctx);
    }

	alphas = FLINT_ARRAY_ALLOC(n, fq_zech_struct);
    Aevals = FLINT_ARRAY_ALLOC(2*(n + 1), fq_zech_mpoly_struct);
    Bevals = Aevals + (n + 1);
	for (i = 0; i < n; i++)
    {
        fq_zech_init(alphas + i, ctx->fqctx);
		fq_zech_mpoly_init(Aevals + i, ctx);
		fq_zech_mpoly_init(Bevals + i, ctx);
    }

    fq_zech_mpoly_init(t1, ctx);
    fq_zech_mpoly_init(t2, ctx);
    fq_zech_mpoly_init(g, ctx);
    fq_zech_mpoly_init(abar, ctx);
    fq_zech_mpoly_init(bbar, ctx);
    fq_zech_mpoly_init(hbar, ctx);

    fq_zech_mpoly_init(Hfac + 0, ctx);
    fq_zech_mpoly_init(Hfac + 1, ctx);
    fq_zech_mpoly_init(Htfac + 0, ctx);
    fq_zech_mpoly_init(Htfac + 1, ctx);

    nmod_mpoly_init(t, smctx);
    nmod_mpoly_init(smHfac + 0, smctx);
    nmod_mpoly_init(smHfac + 1, smctx);

    /* init done */

    _fq_zech_mpoly_set_nmod_mpoly(A, ctx, smA, smctx);
    _fq_zech_mpoly_set_nmod_mpoly(B, ctx, smB, smctx);

    alphas_tries_remaining = 20;

increase_degree:

    edeg++;
    if (edeg > max_degree)
	{
		success = 0;
        goto cleanup;
	}

    fq_zech_mpoly_ctx_change_modulus(ctx, edeg);

    _fq_zech_mpoly_set_nmod_mpoly(A, ctx, smA, smctx);
    _fq_zech_mpoly_set_nmod_mpoly(B, ctx, smB, smctx);

next_alpha:

    if (--alphas_tries_remaining < 0)
	{
		success = 0;
        goto cleanup;
	}

    for (i = 0; i < n; i++)
        fq_zech_rand(alphas + i, state, ctx->fqctx);

    /* ensure deg_X do not drop under evaluation */
    Adegx = fq_zech_mpoly_degree_si(A, 0, ctx);
    Bdegx = fq_zech_mpoly_degree_si(B, 0, ctx);

	for (i = n - 1; i >= 0; i--)
	{
		fq_zech_mpoly_evaluate_one_fq_zech(Aevals + i, i == n - 1 ? A :
                                       Aevals + i + 1, i + 1, alphas + i, ctx);
		fq_zech_mpoly_evaluate_one_fq_zech(Bevals + i, i == n - 1 ? B :
                                       Bevals + i + 1, i + 1, alphas + i, ctx);
		if (Adegx != fq_zech_mpoly_degree_si(Aevals + i, 0, ctx) ||
            Bdegx != fq_zech_mpoly_degree_si(Bevals + i, 0, ctx))
        {
    		goto increase_degree;
        }
	}

    /* univariate gcd */
	success = fq_zech_mpoly_gcd_cofactors(g, abar, bbar,
                                          Aevals + 0, Bevals + 0, ctx) &&
              fq_zech_mpoly_gcd(t1, g, abar, ctx) &&
              fq_zech_mpoly_gcd(t2, g, bbar, ctx);
    if (!success)
        goto cleanup;

    gdegx = fq_zech_mpoly_degree_si(g, 0, ctx);

    if (gdegx == 0)
    {
        /* G is trivial */
        nmod_mpoly_set(Abar, smA, smctx);
        nmod_mpoly_set(Bbar, smB, smctx);
        nmod_mpoly_one(G, smctx);
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
        if (nmod_mpoly_divides(Bbar, smB, smA, smctx))
        {
            nmod_mpoly_set(G, smA, smctx);
            nmod_mpoly_one(Abar, smctx);
            success = 1;
            goto cleanup;
        }

        goto next_alpha;
    }
    else if (gdegx == Bdegx)
    {
        if (nmod_mpoly_divides(Abar, smA, smB, smctx))
        {
            nmod_mpoly_set(G, smB, smctx);
            nmod_mpoly_one(Bbar, smctx);
            success = 1;
            goto cleanup;
        }

        goto next_alpha;
    }

    FLINT_ASSERT(0 < gdegx && gdegx < FLINT_MIN(Adegx, Bdegx));

    /* set Hlcs[n], Glcs[n] (gamma), H, and Hevals */
    if (fq_zech_mpoly_is_one(t1, ctx))
    {
        fq_zech_one(mu1, ctx->fqctx);
        fq_zech_zero(mu2, ctx->fqctx);

        fq_zech_mpoly_swap(hbar, abar, ctx);

        fq_zech_mpolyl_lead_coeff(Hlcs + n, A, 1, ctx);
        fq_zech_mpolyl_lead_coeff(t2, B, 1, ctx);
        success = fq_zech_mpoly_gcd(Glcs + n, Hlcs + n, t2, ctx);
        if (!success)
            goto cleanup;

        H = (fq_zech_mpoly_struct *) A;

        gamma_is_one = fq_zech_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                fq_zech_mpoly_swap(Hevals + i, Aevals + i, ctx);
    }
    else if (fq_zech_mpoly_is_one(t2, ctx))
    {
        fq_zech_zero(mu1, ctx->fqctx);
        fq_zech_one(mu2, ctx->fqctx);

        fq_zech_mpoly_swap(hbar, bbar, ctx);

        fq_zech_mpolyl_lead_coeff(Hlcs + n, B, 1, ctx);
        fq_zech_mpolyl_lead_coeff(t2, A, 1, ctx);
        success = fq_zech_mpoly_gcd(Glcs + n, Hlcs + n, t2, ctx);
        if (!success)
            goto cleanup;

        H = (fq_zech_mpoly_struct *) B;

        gamma_is_one = fq_zech_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                fq_zech_mpoly_swap(Hevals + i, Bevals + i, ctx);
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

        fq_zech_one(mu1, ctx->fqctx);
        fq_zech_rand_not_zero(mu2, state, ctx->fqctx);
        fq_zech_mpoly_scalar_addmul_fq_zech(hbar, abar, bbar, mu2, ctx);

        /* make sure the linear combo did not drop degree */
        if (fq_zech_mpoly_degree_si(hbar, 0, ctx) != FLINT_MAX(Adegx, Bdegx) - gdegx)
            goto next_mu;

        /* make sure the linear combo is prime to g */
        success = fq_zech_mpoly_gcd(t1, hbar, g, ctx);
        if (!success)
            goto cleanup;

        if (!fq_zech_mpoly_is_one(t1, ctx))
            goto next_mu;

        fq_zech_mpolyl_lead_coeff(t1, A, 1, ctx);
        fq_zech_mpolyl_lead_coeff(t2, B, 1, ctx);
        success = fq_zech_mpoly_gcd(Glcs + n, t1, t2, ctx);
        if (!success)
            goto cleanup;

        H = Hevals + n;
        fq_zech_mpoly_scalar_addmul_fq_zech(H, A, B, mu2, ctx);
        fq_zech_mpolyl_lead_coeff(Hlcs + n, H, 1, ctx);

        gamma_is_one = fq_zech_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                fq_zech_mpoly_scalar_addmul_fq_zech(Hevals + i, Aevals + i,
                                                         Bevals + i, mu2, ctx);
    }

    if (!gamma_is_one)
    {
        fq_zech_mpoly_mul(Hevals + n, H, Glcs + n, ctx);
        H = Hevals + n;
        for (i = n - 1; i >= 0; i--)
            fq_zech_mpoly_evaluate_one_fq_zech(Hevals + i, Hevals + i + 1,
                                                       i + 1, alphas + i, ctx);
    }

    success = H->bits <= FLINT_BITS ||
              fq_zech_mpoly_repack_bits_inplace(H, FLINT_BITS, ctx);
    if (!success)
        goto cleanup;

    /* the evals should all fit in H->bits */
    for (i = 0; i < n; i++)
        fq_zech_mpoly_repack_bits_inplace(Hevals + i, H->bits, ctx);

    fq_zech_mpoly_degrees_si(Hdegs, H, ctx);

    /* computed evaluated leading coeffs */
    for (i = n - 1; i >= 0; i--)
    {
        fq_zech_mpoly_evaluate_one_fq_zech(Glcs + i, Glcs + i + 1, i + 1,
                                                              alphas + i, ctx);
        fq_zech_mpoly_evaluate_one_fq_zech(Hlcs + i, Hlcs + i + 1, i + 1,
                                                              alphas + i, ctx);
        /* evaluation could have killed gamma */
        if (fq_zech_mpoly_is_zero(Glcs + i, ctx) ||
            fq_zech_mpoly_is_zero(Hlcs + i, ctx))
        {
            goto increase_degree;
        }
    }

    /* make the leading coefficients match Glcs[0], Hlcs[0] */
    FLINT_ASSERT(fq_zech_mpoly_is_fq_zech(Glcs + 0, ctx) && Glcs[0].length == 1);
    FLINT_ASSERT(fq_zech_mpoly_is_fq_zech(Hlcs + 0, ctx) && Hlcs[0].length == 1);

    fq_zech_inv(q, g->coeffs + 0, ctx->fqctx);
    fq_zech_mul(q, q, Glcs[0].coeffs + 0, ctx->fqctx);
    fq_zech_mpoly_scalar_mul_fq_zech(Hfac + 0, g, q, ctx);

    fq_zech_inv(q, hbar->coeffs + 0, ctx->fqctx);
    fq_zech_mul(q, q, Hlcs[0].coeffs + 0, ctx->fqctx);
    fq_zech_mpoly_scalar_mul_fq_zech(Hfac + 1, hbar, q, ctx);

    for (k = 1; k <= n; k++)
    {
        _fq_zech_mpoly_set_lead0(Htfac + 0, Hfac + 0, Glcs + k, ctx);
        _fq_zech_mpoly_set_lead0(Htfac + 1, Hfac + 1, Hlcs + k, ctx);
        success = fq_zech_mpoly_hlift(k, Htfac, 2, alphas,
                                           k < n ? Hevals + k : H, Hdegs, ctx);
        if (!success)
            goto next_alpha;

        fq_zech_mpoly_swap(Hfac + 0, Htfac + 0, ctx);
        fq_zech_mpoly_swap(Hfac + 1, Htfac + 1, ctx);
    }

    fq_zech_mpoly_scalar_mul_fq_zech(Hfac + 1, Hfac + 1, Hfac[0].coeffs + 0, ctx);
    fq_zech_mpoly_make_monic(Hfac + 0, Hfac + 0, ctx);
    if (!_fq_zech_mpoly_get_nmod_mpoly(smHfac + 0, smctx, Hfac + 0, ctx))
    {
        /* the lifted gcd is not in Fp. possible? */
        goto next_alpha;
    }

    success = nmod_mpolyl_content(t, smHfac + 0, 1, smctx);
    if (!success)
        goto cleanup;

    success = nmod_mpoly_divides(G, smHfac + 0, t, smctx);
    FLINT_ASSERT(success);

    if (fq_zech_is_zero(mu2, ctx->fqctx))
    {
        FLINT_ASSERT(fq_zech_is_one(mu1, ctx->fqctx));

        if (!_fq_zech_mpoly_get_nmod_mpoly(smHfac + 1, smctx, Hfac + 1, ctx))
            goto next_alpha;

        /* the division by t1 should succeed, but let's be careful */
        nmod_mpolyl_lead_coeff(t, G, 1, smctx);
        success = nmod_mpoly_divides(Abar, smHfac + 1, t, smctx) &&
                  nmod_mpoly_divides(Bbar, smB, G, smctx);
    }
    else if (fq_zech_is_zero(mu1, ctx->fqctx))
    {
        FLINT_ASSERT(fq_zech_is_one(mu2, ctx->fqctx));

        if (!_fq_zech_mpoly_get_nmod_mpoly(smHfac + 1, smctx, Hfac + 1, ctx))
            goto next_alpha;

        /* ditto */
        nmod_mpolyl_lead_coeff(t, G, 1, smctx);
        success = nmod_mpoly_divides(Bbar, smHfac + 1, t, smctx) &&
                  nmod_mpoly_divides(Abar, smA, G, smctx);
    }
    else
    {
        FLINT_ASSERT(fq_zech_is_one(mu1, ctx->fqctx));

        success = nmod_mpoly_divides(Abar, smA, G, smctx) &&
                  nmod_mpoly_divides(Bbar, smB, G, smctx);
    }

    if (!success)
        goto next_alpha;

    success = 1;

cleanup:

    fq_zech_clear(q, ctx->fqctx);
    fq_zech_clear(mu1, ctx->fqctx);
    fq_zech_clear(mu2, ctx->fqctx);
    fq_zech_mpoly_clear(A, ctx);
    fq_zech_mpoly_clear(B, ctx);

    flint_randclear(state);

    flint_free(Hdegs);

    for (i = 0; i < n + 1; i++)
    {
        fq_zech_mpoly_clear(Glcs + i, ctx);
        fq_zech_mpoly_clear(Hlcs + i, ctx);
        fq_zech_mpoly_clear(Hevals + i, ctx);
    }
    flint_free(Glcs);

	for (i = 0; i < n; i++)
    {
        fq_zech_clear(alphas + i, ctx->fqctx);
		fq_zech_mpoly_clear(Aevals + i, ctx);
		fq_zech_mpoly_clear(Bevals + i, ctx);
    }
    flint_free(alphas);
    flint_free(Aevals);

    fq_zech_mpoly_clear(t1, ctx);
    fq_zech_mpoly_clear(t2, ctx);
    fq_zech_mpoly_clear(g, ctx);
    fq_zech_mpoly_clear(abar, ctx);
    fq_zech_mpoly_clear(bbar, ctx);
    fq_zech_mpoly_clear(hbar, ctx);

    fq_zech_mpoly_clear(Hfac + 0, ctx);
    fq_zech_mpoly_clear(Hfac + 1, ctx);
    fq_zech_mpoly_clear(Htfac + 0, ctx);
    fq_zech_mpoly_clear(Htfac + 1, ctx);

    nmod_mpoly_clear(t, smctx);
    nmod_mpoly_clear(smHfac + 0, smctx);
    nmod_mpoly_clear(smHfac + 1, smctx);

    fq_zech_mpoly_ctx_clear(ctx);

    if (success)
    {
        nmod_mpoly_repack_bits_inplace(G, bits, smctx);
        nmod_mpoly_repack_bits_inplace(Abar, bits, smctx);
        nmod_mpoly_repack_bits_inplace(Bbar, bits, smctx);

        FLINT_ASSERT(G->length > 0);
        FLINT_ASSERT(Abar->length > 0);
        FLINT_ASSERT(Bbar->length > 0);
    }

	return success;
}

/* should find its way back here in interesting cases */
int nmod_mpoly_gcd_hensel(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    if (nmod_mpoly_is_zero(A, ctx) || nmod_mpoly_is_zero(B, ctx))
        return nmod_mpoly_gcd(G, A, B, ctx);

    return _nmod_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_HENSEL);
}

