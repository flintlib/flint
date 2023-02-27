/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


int fmpz_mpoly_factor_irred_wang(
    fmpz_mpolyv_t fac,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_factor_t lcAfac,
    int lcAfac_irred,
    const fmpz_mpoly_t lcA,
    const fmpz_mpoly_ctx_t ctx,
    flint_rand_t state,
    zassenhaus_prune_t zas,
    int allow_shift)
{
    int success, kfails = 0;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, r;
    fmpz * alpha;
    slong alpha_modulus, alpha_count;
    fmpz_mpoly_struct * Aevals;
    slong * degs, * tdegs;
    fmpz_mpolyv_t tfac;
    fmpz_mpoly_t t, Acopy;
    fmpz_mpoly_struct * newA;
    fmpz_poly_t Au;
    fmpz_poly_factor_t Aufac;
    fmpz_mpoly_t m, mpow;
    fmpz_mpolyv_t new_lcs, lc_divs;
    fmpz_t q;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(fmpz_mpoly_factor_matches(lcA, lcAfac, ctx));

    fmpz_init(q);

    fmpz_mpoly_init(Acopy, ctx);
    fmpz_mpoly_init(m, ctx);
    fmpz_mpoly_init(mpow, ctx);

    fmpz_mpolyv_init(new_lcs, ctx);
    fmpz_mpolyv_init(lc_divs, ctx);

    fmpz_poly_factor_init(Aufac);
    fmpz_poly_init(Au);

    degs  = (slong *) flint_malloc(2*(n + 1)*sizeof(slong));
    tdegs = degs + (n + 1);
    alpha = _fmpz_vec_init(n);
    Aevals    = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < n; i++)
        fmpz_mpoly_init(Aevals + i, ctx);

    fmpz_mpolyv_init(tfac, ctx);
    fmpz_mpoly_init(t, ctx);

    /* init done */

    fmpz_mpoly_degrees_si(degs, A, ctx);

    alpha_count = 0;
    alpha_modulus = 1;
    goto got_alpha;

next_alpha:

    if (!allow_shift)
    {
        success = 0;
        goto cleanup;
    }

    alpha_count++;
    if (alpha_count >= alpha_modulus)
    {
        alpha_count = 0;
        alpha_modulus++;
        if (alpha_modulus/1024 > ctx->minfo->nvars)
        {
            success = 0;
            goto cleanup;
        }
    }

    for (i = 0; i < n; i++)
        fmpz_set_si(alpha + i, n_urandint(state, alpha_modulus) - alpha_modulus/2);

got_alpha:

#if FLINT_WANT_ASSERT
    fmpz_mpoly_degrees_si(tdegs, A, ctx);
    for (i = 0; i < n + 1; i++)
        FLINT_ASSERT(degs[i] == tdegs[i]);
#endif

    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
        fmpz_mpoly_degrees_si(tdegs, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
        {
            if (tdegs[j] != degs[j])
                goto next_alpha;
        }
    }

    /* make sure our univar is squarefree */
    FLINT_ASSERT(fmpz_mpoly_is_fmpz_poly(Aevals + 0, 0, ctx));
    success = fmpz_mpoly_get_fmpz_poly(Au, Aevals + 0, 0, ctx);
    FLINT_ASSERT(success);
    fmpz_poly_factor(Aufac, Au);
    r = Aufac->num;

    zassenhaus_prune_start_add_factors(zas);
    for (j = 0; j < r; j++)
        zassenhaus_prune_add_factor(zas, fmpz_poly_degree(Aufac->p + j),
                                                                Aufac->exp[j]);
    zassenhaus_prune_end_add_factors(zas);

    if ((r < 2 && Aufac->exp[0] == 1) ||
        zassenhaus_prune_must_be_irreducible(zas))
    {
        fmpz_mpolyv_fit_length(fac, 1, ctx);
        fac->length = 1;
        fmpz_mpoly_set(fac->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    for (j = 0; j < r; j++)
    {
        if (Aufac->exp[j] != 1)
            goto next_alpha;
    }

    fmpz_mpolyv_fit_length(lc_divs, r, ctx);
    lc_divs->length = r;
    if (lcAfac->num > 0)
    {
        success = 0;
        if (lcAfac_irred)
            success = fmpz_mpoly_factor_lcc_wang(lc_divs->coeffs, lcAfac,
                                           &Aufac->c, Aufac->p, r, alpha, ctx);
        if (!success)
        {
            success = fmpz_mpoly_factor_lcc_kaltofen(lc_divs->coeffs, lcAfac,
                                                A, r, alpha, degs, Aufac, ctx);
            if (success < 0 || (success == 0 && ++kfails < 4))
                goto next_alpha;
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            FLINT_ASSERT(Aufac->p[i].length > 0);
            fmpz_mpoly_set_fmpz(lc_divs->coeffs + i,
                             Aufac->p[i].coeffs + Aufac->p[i].length - 1, ctx);
        }
    }

    FLINT_ASSERT(r > 1);
    success = fmpz_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    if (!success)
        goto next_alpha;
    for (i = 1; i < r; i++)
    {
        success = fmpz_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        if (!success)
            goto next_alpha;
    }

    fmpz_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fmpz_mpoly_is_one(mpow, ctx))
    {
        newA = (fmpz_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fmpz_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = -1;
        goto cleanup;
    }

    fmpz_mpoly_degrees_si(tdegs, newA, ctx);

    fmpz_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(t, mpow, i + 1, alpha + i, ctx);
        fmpz_mpoly_swap(t, mpow, ctx);
        fmpz_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
    }

    fmpz_mpolyv_fit_length(new_lcs, (n + 1)*r, ctx);
    i = n;
    for (j = 0; j < r; j++)
    {
        fmpz_mpoly_mul(new_lcs->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_evaluate_one_fmpz(new_lcs->coeffs + i*r + j,
                       new_lcs->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fmpz_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(fmpz_mpoly_is_fmpz(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fmpz_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        FLINT_ASSERT(fmpz_divisible(new_lcs->coeffs[i].coeffs + 0, Aufac->p[i].coeffs + Aufac->p[i].length - 1));

        fmpz_divexact(q, new_lcs->coeffs[i].coeffs + 0,
                                  Aufac->p[i].coeffs + Aufac->p[i].length - 1);
        _fmpz_mpoly_set_fmpz_poly(fac->coeffs + i, newA->bits,
                               Aufac->p[i].coeffs, Aufac->p[i].length, 0, ctx);
        fmpz_mpoly_scalar_mul_fmpz(fac->coeffs + i, fac->coeffs + i, q, ctx);
    }

    fmpz_mpolyv_fit_length(tfac, r, ctx);
    tfac->length = r;
    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            _fmpz_mpoly_set_lead0(tfac->coeffs + i, fac->coeffs + i,
                                               new_lcs->coeffs + k*r + i, ctx);
        }

        success = fmpz_mpoly_hlift(k, tfac->coeffs, r, alpha,
                                      k < n ? Aevals + k : newA, tdegs, ctx);
        if (!success)
            goto next_alpha;

        fmpz_mpolyv_swap(tfac, fac, ctx);
    }

    if (fmpz_mpoly_is_fmpz(m, ctx))
    {
        for (i = 0; i < r; i++)
        {
            _fmpz_vec_content(q, fac->coeffs[i].coeffs, fac->coeffs[i].length);
            if (fmpz_sgn(fac->coeffs[i].coeffs + 0) < 0)
                fmpz_neg(q, q);
            fmpz_mpoly_scalar_divexact_fmpz(fac->coeffs + i,
                                            fac->coeffs + i, q, ctx);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            /* hlift should not have returned any large bits */
            FLINT_ASSERT(fac->coeffs[i].bits <= FLINT_BITS);

            if (!fmpz_mpolyl_content(t, fac->coeffs + i, 1, ctx))
            {
                success = -1;
                goto cleanup;
            }

            success = fmpz_mpoly_divides(fac->coeffs + i, fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
            fmpz_mpoly_unit_normalize(fac->coeffs + i, ctx);
        }
    }

    success = 1;

cleanup:

    fmpz_clear(q);

    fmpz_mpolyv_clear(new_lcs, ctx);
    fmpz_mpolyv_clear(lc_divs, ctx);

    fmpz_poly_factor_clear(Aufac);
        
    _fmpz_vec_clear(alpha, n);

    for (i = 0; i < n; i++)
        fmpz_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);

    flint_free(degs); /* and tdegs */
    fmpz_mpolyv_clear(tfac, ctx);
    fmpz_mpoly_clear(t, ctx);

    fmpz_mpoly_clear(Acopy, ctx);
    fmpz_mpoly_clear(m, ctx);
    fmpz_mpoly_clear(mpow, ctx);

    fmpz_poly_clear(Au);

#if FLINT_WANT_ASSERT
    if (success)
    {
        fmpz_mpoly_t prod;
        fmpz_mpoly_init(prod, ctx);
        fmpz_mpoly_one(prod, ctx);
        for (i = 0; i < fac->length; i++)
            fmpz_mpoly_mul(prod, prod, fac->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mpoly_equal(prod, A, ctx));
        fmpz_mpoly_clear(prod, ctx);
    }
#endif

    return success;
}

