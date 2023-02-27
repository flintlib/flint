/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

/*
    return 1: success
           0: failed
          -1: exception (large exps)
*/
int nmod_mpoly_factor_irred_smprime_wang(
    nmod_mpolyv_t fac,
    const nmod_mpoly_t A,
    const nmod_mpoly_factor_t lcAfac,
    const nmod_mpoly_t lcA,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    int alphas_tries_remaining, alphabetas_tries_remaining, alphabetas_length;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, r;
    mp_limb_t * alpha;
    n_poly_struct * alphabetas;
    nmod_mpoly_struct * Aevals;
    slong * degs, * degeval;
    nmod_mpolyv_t tfac;
    nmod_mpoly_t t, Acopy;
    nmod_mpoly_struct * newA;
    n_poly_t Abfc;
    n_bpoly_t Ab;
    n_tpoly_t Abfp;
    nmod_mpoly_t m, mpow;
    nmod_mpolyv_t new_lcs, lc_divs;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 1);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    nmod_mpoly_init(Acopy, ctx);
    nmod_mpoly_init(m, ctx);
    nmod_mpoly_init(mpow, ctx);

    nmod_mpolyv_init(new_lcs, ctx);
    nmod_mpolyv_init(lc_divs, ctx);

    n_poly_init(Abfc);
    n_tpoly_init(Abfp);
    n_bpoly_init(Ab);

    degs    = (slong *) flint_malloc((n + 1)*sizeof(slong));
    degeval = (slong *) flint_malloc((n + 1)*sizeof(slong));
    alpha   = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));
    alphabetas = (n_poly_struct *) flint_malloc(n*sizeof(n_poly_struct));
    Aevals  = (nmod_mpoly_struct *) flint_malloc(n*sizeof(nmod_mpoly_struct));
    for (i = 0; i < n; i++)
    {
        n_poly_init(alphabetas + i);
        nmod_mpoly_init(Aevals + i, ctx);
    }
    nmod_mpolyv_init(tfac, ctx);
    nmod_mpoly_init(t, ctx);

    /* init done */

    alphabetas_length = 2;
    alphas_tries_remaining = 10;
    nmod_mpoly_degrees_si(degs, A, ctx);

next_alpha:

    if (--alphas_tries_remaining < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < n; i++)
        alpha[i] = n_urandint(state, ctx->mod.n - 1) + 1;

    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        nmod_mpoly_evaluate_one_ui(Aevals + i,
                        i == n - 1 ? A : Aevals + i + 1, i + 1, alpha[i], ctx);
        nmod_mpoly_degrees_si(degeval, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
            if (degeval[j] != degs[j])
                goto next_alpha;
    }

    /* make sure univar is squarefree */
    nmod_mpoly_derivative(t, Aevals + 0, 0, ctx);
    if (!nmod_mpoly_gcd(t, t, Aevals + 0, ctx))
    {
        success = -1;
        goto cleanup;
    }
    if (!nmod_mpoly_is_one(t, ctx))
        goto next_alpha;

    alphabetas_tries_remaining = 2 + alphabetas_length;

next_alphabetas:

    if (--alphabetas_tries_remaining < 0)
    {
        if (++alphabetas_length > 10)
        {
            success = 0;
            goto cleanup;
        }
        goto next_alpha;
    }

    for (i = 0; i < n; i++)
    {
        n_poly_fit_length(alphabetas + i, alphabetas_length);
        alphabetas[i].coeffs[0] = alpha[i];
        for (j = 1; j < alphabetas_length; j++)
            alphabetas[i].coeffs[j] = n_urandint(state, ctx->mod.n);
        alphabetas[i].length = alphabetas_length;
        _n_poly_normalise(alphabetas + i);
    }

    _nmod_mpoly_eval_rest_to_n_bpoly(Ab, A, alphabetas, ctx);
    success = n_bpoly_mod_factor_smprime(Abfc, Abfp, Ab, 0, ctx->mod);
    if (!success)
    {
        FLINT_ASSERT(0 && "this should not happen");
        goto next_alpha;
    }

    r = Abfp->length;

    if (r < 2)
    {
        nmod_mpolyv_fit_length(fac, 1, ctx);
        fac->length = 1;
        nmod_mpoly_set(fac->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    nmod_mpolyv_fit_length(lc_divs, r, ctx);
    lc_divs->length = r;
    if (lcAfac->num > 0)
    {
        success = nmod_mpoly_factor_lcc_wang(lc_divs->coeffs, lcAfac,
                                       Abfc, Abfp->coeffs, r, alphabetas, ctx);
        if (!success)
            goto next_alphabetas;
    }
    else
    {
        for (i = 0; i < r; i++)
            nmod_mpoly_one(lc_divs->coeffs + i, ctx);
    }

    success = nmod_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < r; i++)
    {
        success = nmod_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }

    nmod_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (nmod_mpoly_is_one(mpow, ctx))
    {
        newA = (nmod_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        nmod_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    nmod_mpoly_degrees_si(degs, newA, ctx);

    nmod_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        nmod_mpoly_evaluate_one_ui(t, mpow, i + 1, alpha[i], ctx);
        nmod_mpoly_swap(t, mpow, ctx);
        nmod_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
        nmod_mpoly_repack_bits_inplace(Aevals + i, newA->bits, ctx);
    }

    nmod_mpolyv_fit_length(new_lcs, (n + 1)*r, ctx);
    i = n;
    for (j = 0; j < r; j++)
    {
        nmod_mpoly_mul(new_lcs->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            nmod_mpoly_evaluate_one_ui(new_lcs->coeffs + i*r + j,
                        new_lcs->coeffs + (i + 1)*r + j, i + 1, alpha[i], ctx);
        }
    }

    nmod_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        mp_limb_t q;
        FLINT_ASSERT(nmod_mpoly_is_ui(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(nmod_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        _nmod_mpoly_set_n_bpoly_var1_zero(fac->coeffs + i, newA->bits,
                                                     Abfp->coeffs + i, 0, ctx);
        FLINT_ASSERT(fac->coeffs[i].length > 0);
        q = nmod_inv(fac->coeffs[i].coeffs[0], ctx->mod);
        q = nmod_mul(q, new_lcs->coeffs[0*r + i].coeffs[0], ctx->mod);
        nmod_mpoly_scalar_mul_nmod_invertible(fac->coeffs + i,
                                                      fac->coeffs + i, q, ctx);
    }

    nmod_mpolyv_fit_length(tfac, r, ctx);
    tfac->length = r;
    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            _nmod_mpoly_set_lead0(tfac->coeffs + i, fac->coeffs + i,
                                               new_lcs->coeffs + k*r + i, ctx);
        }

        success = nmod_mpoly_hlift(k, tfac->coeffs, r, alpha,
                                         k < n ? Aevals + k : newA, degs, ctx);
        if (!success)
            goto next_alphabetas;

        nmod_mpolyv_swap(tfac, fac, ctx);
    }

    if (!nmod_mpoly_is_ui(m, ctx))
    {
        for (i = 0; i < r; i++)
        {
            /* hlift should not have returned any large bits */
            FLINT_ASSERT(fac->coeffs[i].bits <= FLINT_BITS);

            if (!nmod_mpolyl_content(t, fac->coeffs + i, 1, ctx))
            {
                success = -1;
                goto cleanup;
            }

            success = nmod_mpoly_divides(fac->coeffs + i, fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
        }
    }

    for (i = 0; i < r; i++)
    {
        /* hlift should not have returned any large bits */
        FLINT_ASSERT(fac->coeffs[i].bits <= FLINT_BITS);
        nmod_mpoly_make_monic(fac->coeffs + i, fac->coeffs + i, ctx);
    }

    success = 1;

cleanup:

    nmod_mpolyv_clear(new_lcs, ctx);
    nmod_mpolyv_clear(lc_divs, ctx);

    n_poly_clear(Abfc);
    n_tpoly_clear(Abfp);
    n_bpoly_clear(Ab);

    for (i = 0; i < n; i++)
    {
        nmod_mpoly_clear(Aevals + i, ctx);
        n_poly_clear(alphabetas + i);
    }
    flint_free(alphabetas);
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(degs);
    flint_free(degeval);

    nmod_mpolyv_clear(tfac, ctx);
    nmod_mpoly_clear(t, ctx);

    nmod_mpoly_clear(Acopy, ctx);
    nmod_mpoly_clear(m, ctx);
    nmod_mpoly_clear(mpow, ctx);

#if FLINT_WANT_ASSERT
    if (success)
    {
        nmod_mpoly_t prod;
        nmod_mpoly_init(prod, ctx);
        nmod_mpoly_one(prod, ctx);
        for (i = 0; i < fac->length; i++)
            nmod_mpoly_mul(prod, prod, fac->coeffs + i, ctx);
        FLINT_ASSERT(nmod_mpoly_equal(prod, A, ctx));
        nmod_mpoly_clear(prod, ctx);
    }
#endif

    return success;
}
