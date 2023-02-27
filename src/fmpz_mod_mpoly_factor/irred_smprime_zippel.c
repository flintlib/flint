/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"

/*
    return 1: success
           0: failed
          -1: exception (large exps)
*/
int fmpz_mod_mpoly_factor_irred_smprime_zippel(
    fmpz_mod_mpolyv_t fac,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_factor_t lcAfac,
    const fmpz_mod_mpoly_t lcA,
    const fmpz_mod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    int alphas_tries_remaining, alphabetas_tries_remaining, alphabetas_length;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, r;
    fmpz * alpha;
    fmpz_mod_poly_struct * alphabetas;
    fmpz_mod_mpoly_struct * Aevals;
    slong * degs, * degeval;
    fmpz_mod_mpolyv_t tfac;
    fmpz_mod_mpoly_t t, Acopy;
    fmpz_mod_mpoly_struct * newA;
    fmpz_mod_poly_t Abfc;
    fmpz_mod_bpoly_t Ab;
    fmpz_mod_tpoly_t Abfp;
    fmpz_mod_mpoly_t m, mpow;
    fmpz_mod_mpolyv_t new_lcs, lc_divs;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 1);
    FLINT_ASSERT(fmpz_is_one(A->coeffs + 0));
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    fmpz_mod_mpoly_init(Acopy, ctx);
    fmpz_mod_mpoly_init(m, ctx);
    fmpz_mod_mpoly_init(mpow, ctx);

    fmpz_mod_mpolyv_init(new_lcs, ctx);
    fmpz_mod_mpolyv_init(lc_divs, ctx);

    fmpz_mod_poly_init(Abfc, ctx->ffinfo);
    fmpz_mod_tpoly_init(Abfp, ctx->ffinfo);
    fmpz_mod_bpoly_init(Ab, ctx->ffinfo);

    degs    = FLINT_ARRAY_ALLOC(n + 1, slong);
    degeval = FLINT_ARRAY_ALLOC(n + 1, slong);
    alpha   = _fmpz_vec_init(n);
    alphabetas = FLINT_ARRAY_ALLOC(n, fmpz_mod_poly_struct);
    Aevals  = FLINT_ARRAY_ALLOC(n, fmpz_mod_mpoly_struct);
    for (i = 0; i < n; i++)
    {
        fmpz_mod_poly_init(alphabetas + i, ctx->ffinfo);
        fmpz_mod_mpoly_init(Aevals + i, ctx);
    }
    fmpz_mod_mpolyv_init(tfac, ctx);
    fmpz_mod_mpoly_init(t, ctx);

    /* init done */

    alphabetas_length = 2;
    alphas_tries_remaining = 10;
    fmpz_mod_mpoly_degrees_si(degs, A, ctx);

next_alpha:

    if (--alphas_tries_remaining < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < n; i++)
        fmpz_mod_rand_not_zero(alpha + i, state, ctx->ffinfo);

    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mod_mpoly_evaluate_one_fmpz(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
        fmpz_mod_mpoly_degrees_si(degeval, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
            if (degeval[j] != degs[j])
                goto next_alpha;
    }

    /* make sure univar is squarefree */
    fmpz_mod_mpoly_derivative(t, Aevals + 0, 0, ctx);
    if (!fmpz_mod_mpoly_gcd(t, t, Aevals + 0, ctx))
    {
        success = -1;
        goto cleanup;
    }
    if (!fmpz_mod_mpoly_is_one(t, ctx))
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
        fmpz_mod_poly_fit_length(alphabetas + i, alphabetas_length, ctx->ffinfo);
        fmpz_set(alphabetas[i].coeffs + 0, alpha + i);
        for (j = 1; j < alphabetas_length; j++)
            fmpz_mod_rand(alphabetas[i].coeffs + j, state, ctx->ffinfo);
        _fmpz_mod_poly_set_length(alphabetas + i, alphabetas_length);
        _fmpz_mod_poly_normalise(alphabetas + i);
    }

    _fmpz_mod_mpoly_eval_rest_to_fmpz_mod_bpoly(Ab, A, alphabetas, ctx);
    success = fmpz_mod_bpoly_factor_smprime(Abfc, Abfp, Ab, 0, ctx->ffinfo);
    if (!success)
    {
        FLINT_ASSERT(0 && "this should not happen");
        goto next_alpha;
    }

    r = Abfp->length;

    if (r < 2)
    {
        fmpz_mod_mpolyv_fit_length(fac, 1, ctx);
        fac->length = 1;
        fmpz_mod_mpoly_set(fac->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    fmpz_mod_mpolyv_fit_length(lc_divs, r, ctx);
    lc_divs->length = r;
    if (lcAfac->num > 0)
    {
        success = fmpz_mod_mpoly_factor_lcc_wang(lc_divs->coeffs, lcAfac,
                                       Abfc, Abfp->coeffs, r, alphabetas, ctx);
        if (!success)
            goto next_alphabetas;
    }
    else
    {
        for (i = 0; i < r; i++)
            fmpz_mod_mpoly_one(lc_divs->coeffs + i, ctx);
    }

    success = fmpz_mod_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < r; i++)
    {
        success = fmpz_mod_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }

    fmpz_mod_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fmpz_mod_mpoly_is_one(mpow, ctx))
    {
        newA = (fmpz_mod_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fmpz_mod_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    fmpz_mod_mpoly_degrees_si(degs, newA, ctx);

    for (i = 0; i < n + 1; i++)
    {
        if (FLINT_BIT_COUNT(degs[i]) >= FLINT_BITS/3)
        {
            success = -1;
            goto cleanup;
        }
    }

    fmpz_mod_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mod_mpoly_evaluate_one_fmpz(t, mpow, i + 1, alpha + i, ctx);
        fmpz_mod_mpoly_swap(t, mpow, ctx);
        fmpz_mod_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
        fmpz_mod_mpoly_repack_bits_inplace(Aevals + i, newA->bits, ctx);
    }

    fmpz_mod_mpolyv_fit_length(new_lcs, (n + 1)*r, ctx);
    i = n;
    for (j = 0; j < r; j++)
    {
        fmpz_mod_mpoly_mul(new_lcs->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fmpz_mod_mpoly_evaluate_one_fmpz(new_lcs->coeffs + i*r + j,
                       new_lcs->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fmpz_mod_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        fmpz_t q;
        fmpz_init(q);
        FLINT_ASSERT(fmpz_mod_mpoly_is_fmpz(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fmpz_mod_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        _fmpz_mod_mpoly_set_fmpz_mod_bpoly_var1_zero(fac->coeffs + i, newA->bits,
                                                     Abfp->coeffs + i, 0, ctx);
        FLINT_ASSERT(fac->coeffs[i].length > 0);
        fmpz_mod_inv(q, fac->coeffs[i].coeffs + 0, ctx->ffinfo);
        fmpz_mod_mul(q, q, new_lcs->coeffs[0*r + i].coeffs + 0, ctx->ffinfo);
        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fac->coeffs + i,
                                                      fac->coeffs + i, q, ctx);
        fmpz_clear(q);
    }

    fmpz_mod_mpolyv_fit_length(tfac, r, ctx);
    tfac->length = r;
    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            _fmpz_mod_mpoly_set_lead0(tfac->coeffs + i, fac->coeffs + i,
                                               new_lcs->coeffs + k*r + i, ctx);
        }

        if (k > 2)
        {
            success = fmpz_mod_mpoly_hlift_zippel(k, tfac->coeffs, r, alpha,
                                  k < n ? Aevals + k : newA, degs, ctx, state);
        }
        else
        {
            success = fmpz_mod_mpoly_hlift(k, tfac->coeffs, r, alpha,
                                         k < n ? Aevals + k : newA, degs, ctx);
        }

        if (!success)
            goto next_alphabetas;

        fmpz_mod_mpolyv_swap(tfac, fac, ctx);
    }

    if (!fmpz_mod_mpoly_is_fmpz(m, ctx))
    {
        for (i = 0; i < r; i++)
        {
            /* hlift should not have returned any large bits */
            FLINT_ASSERT(fac->coeffs[i].bits <= FLINT_BITS);

            if (!fmpz_mod_mpolyl_content(t, fac->coeffs + i, 1, ctx))
            {
                success = -1;
                goto cleanup;
            }

            success = fmpz_mod_mpoly_divides(fac->coeffs + i, fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
        }
    }

    for (i = 0; i < r; i++)
    {
        /* hlift should not have returned any large bits */
        FLINT_ASSERT(fac->coeffs[i].bits <= FLINT_BITS);
        fmpz_mod_mpoly_make_monic(fac->coeffs + i, fac->coeffs + i, ctx);
    }

    success = 1;

cleanup:

    fmpz_mod_mpolyv_clear(new_lcs, ctx);
    fmpz_mod_mpolyv_clear(lc_divs, ctx);

    fmpz_mod_poly_clear(Abfc, ctx->ffinfo);
    fmpz_mod_tpoly_clear(Abfp, ctx->ffinfo);
    fmpz_mod_bpoly_clear(Ab, ctx->ffinfo);

    for (i = 0; i < n; i++)
    {
        fmpz_mod_mpoly_clear(Aevals + i, ctx);
        fmpz_mod_poly_clear(alphabetas + i, ctx->ffinfo);
    }
    flint_free(alphabetas);
    flint_free(Aevals);
    _fmpz_vec_clear(alpha, n);
    flint_free(degs);
    flint_free(degeval);

    fmpz_mod_mpolyv_clear(tfac, ctx);
    fmpz_mod_mpoly_clear(t, ctx);

    fmpz_mod_mpoly_clear(Acopy, ctx);
    fmpz_mod_mpoly_clear(m, ctx);
    fmpz_mod_mpoly_clear(mpow, ctx);

#if FLINT_WANT_ASSERT
    if (success)
    {
        fmpz_mod_mpoly_t prod;
        fmpz_mod_mpoly_init(prod, ctx);
        fmpz_mod_mpoly_one(prod, ctx);
        for (i = 0; i < fac->length; i++)
            fmpz_mod_mpoly_mul(prod, prod, fac->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mod_mpoly_equal(prod, A, ctx));
        fmpz_mod_mpoly_clear(prod, ctx);
    }
#endif

FLINT_ASSERT(success);
    return success;
}
