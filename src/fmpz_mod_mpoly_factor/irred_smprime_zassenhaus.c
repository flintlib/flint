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
    return:
        1: success
        0: lift is impossible
       -1: failed, don't try again
*/
static int _try_lift(
    fmpz_mod_mpolyv_t qfac,
    const fmpz_mod_mpoly_t q,
    const fmpz_mod_mpolyv_t pfac,
    const fmpz_mod_mpoly_t p,
    slong m,
    const fmpz * alphas,
    slong n,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    slong * newdeg;
    fmpz_mod_mpoly_t lcq, lcp, t, newq;

    FLINT_ASSERT(pfac->length > 1);

    newdeg = FLINT_ARRAY_ALLOC(n + 1, slong);
    fmpz_mod_mpoly_init(lcq, ctx);
    fmpz_mod_mpoly_init(lcp, ctx);
    fmpz_mod_mpoly_init(t, ctx);
    fmpz_mod_mpoly_init(newq, ctx);

#if FLINT_WANT_ASSERT
    fmpz_mod_mpoly_one(t, ctx);
    for (i = 0; i < pfac->length; i++)
        fmpz_mod_mpoly_mul(t, t, pfac->coeffs + i, ctx);
    FLINT_ASSERT(fmpz_mod_mpoly_equal(t, p, ctx));
#endif

    _fmpz_mod_mpoly_get_lead0(lcq, q, ctx);
    fmpz_mod_mpoly_evaluate_one_fmpz(lcp, lcq, m, alphas + m - 1, ctx);

    FLINT_ASSERT(lcp->length > 0);

    fmpz_mod_mpoly_pow_ui(t, lcq, pfac->length - 1, ctx);
    fmpz_mod_mpoly_mul(newq, q, t, ctx);

    if (newq->bits > FLINT_BITS)
    {
        success = -1;
        goto cleanup;
    }

    fmpz_mod_mpoly_degrees_si(newdeg, newq, ctx);

    fmpz_mod_mpolyv_fit_length(qfac, pfac->length, ctx);
    qfac->length = pfac->length;
    for (i = 0; i < pfac->length; i++)
    {
        _fmpz_mod_mpoly_get_lead0(t, pfac->coeffs + i, ctx);
        success = fmpz_mod_mpoly_divides(t, lcp, t, ctx);
        FLINT_ASSERT(success);
        fmpz_mod_mpoly_mul(qfac->coeffs + i, pfac->coeffs + i, t, ctx);
        _fmpz_mod_mpoly_set_lead0(qfac->coeffs + i, qfac->coeffs + i, lcq, ctx);
    }

    success = fmpz_mod_mpoly_hlift(m, qfac->coeffs, qfac->length,
                                                    alphas, newq, newdeg, ctx);
    if (!success)
        goto cleanup;

    for (i = 0; i < qfac->length; i++)
    {
        /* hlift should not have returned any large bits */
        FLINT_ASSERT(qfac->coeffs[i].bits <= FLINT_BITS);

        if (!fmpz_mod_mpolyl_content(t, qfac->coeffs + i, 1, ctx))
        {
            success = -1;
            goto cleanup;
        }

        success = fmpz_mod_mpoly_divides(qfac->coeffs + i, qfac->coeffs + i, t, ctx);
        FLINT_ASSERT(success);
        fmpz_mod_mpoly_make_monic(qfac->coeffs + i, qfac->coeffs + i, ctx);
    }

    success = 1;

cleanup:

    flint_free(newdeg);
    fmpz_mod_mpoly_clear(lcq, ctx);
    fmpz_mod_mpoly_clear(lcp, ctx);
    fmpz_mod_mpoly_clear(t, ctx);
    fmpz_mod_mpoly_clear(newq, ctx);

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_one(t, ctx);
        for (i = 0; i < qfac->length; i++)
            fmpz_mod_mpoly_mul(t, t, qfac->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mod_mpoly_equal(t, q, ctx));
        fmpz_mod_mpoly_clear(t, ctx);
    }
#endif

    return success;
}


int fmpz_mod_mpoly_factor_irred_smprime_zassenhaus(
    fmpz_mod_mpolyv_t fac,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    int tries_remaining = 10;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, m, len;
    slong * subset;
    fmpz * alphas;
    fmpz_mod_mpoly_struct * Aevals;
    slong * deg, * degeval;
    fmpz_mod_mpolyv_t qfac, pfac, tfac, dfac;
    fmpz_mod_mpoly_t t, p, q;
    fmpz_mod_poly_t c;
    fmpz_mod_bpoly_t B;
    fmpz_mod_tpoly_t F;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 1);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    subset = FLINT_ARRAY_ALLOC(4, slong);
    alphas = _fmpz_vec_init(n);
    Aevals = FLINT_ARRAY_ALLOC(n, fmpz_mod_mpoly_struct);
    deg     = FLINT_ARRAY_ALLOC(n + 1, slong);
    degeval = FLINT_ARRAY_ALLOC(n + 1, slong);
    for (i = 0; i < n; i++)
        fmpz_mod_mpoly_init(Aevals + i, ctx);
    fmpz_mod_mpolyv_init(pfac, ctx);
    fmpz_mod_mpolyv_init(qfac, ctx);
    fmpz_mod_mpolyv_init(tfac, ctx);
    fmpz_mod_mpolyv_init(dfac, ctx);
    fmpz_mod_mpoly_init(t, ctx);
    fmpz_mod_mpoly_init(p, ctx);
    fmpz_mod_mpoly_init(q, ctx);
    fmpz_mod_poly_init(c, ctx->ffinfo);
    fmpz_mod_bpoly_init(B, ctx->ffinfo);
    fmpz_mod_tpoly_init(F, ctx->ffinfo);

    fmpz_mod_mpoly_degrees_si(deg, A, ctx);

next_alpha:

    if (--tries_remaining < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < n; i++)
        fmpz_mod_rand_not_zero(alphas + i, state, ctx->ffinfo);

    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mod_mpoly_evaluate_one_fmpz(Aevals + i, i == n - 1 ? A :
                                       Aevals + i + 1, i + 1, alphas + i, ctx);
        fmpz_mod_mpoly_repack_bits_inplace(Aevals + i, A->bits, ctx);
        fmpz_mod_mpoly_degrees_si(degeval, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
            if (degeval[j] != deg[j])
                goto next_alpha;
    }

    /* make sure univar is squarefree */
    fmpz_mod_mpoly_derivative(t, Aevals + 0, 0, ctx);
    success = fmpz_mod_mpoly_gcd(t, t, Aevals + 0, ctx);
    if (!success)
        goto cleanup;
    if (!fmpz_mod_mpoly_is_one(t, ctx))
        goto next_alpha;

    /* make evaluations primitive */
    for (i = n - 1; i > 0; i--)
    {
        success = fmpz_mod_mpolyl_content(t, Aevals + i, 1, ctx);
        if (!success)
            goto cleanup;
        success = fmpz_mod_mpoly_divides(Aevals + i, Aevals + i, t, ctx);
        FLINT_ASSERT(success);
        fmpz_mod_mpoly_make_monic(Aevals + i, Aevals + i, ctx);
        fmpz_mod_mpoly_repack_bits_inplace(Aevals + i, A->bits, ctx);
    }

    fmpz_mod_mpoly_get_fmpz_mod_bpoly(B, Aevals + 1, 0, 1, ctx);
    success = fmpz_mod_bpoly_factor_smprime(c, F, B, 1, ctx->ffinfo);
    if (!success)
        goto next_alpha;

    FLINT_ASSERT(fmpz_mod_poly_degree(c, ctx->ffinfo) == 0);

    fmpz_mod_mpolyv_fit_length(pfac, F->length, ctx);
    pfac->length = F->length;
    for (i = 0; i < F->length; i++)
    {
        fmpz_mod_mpoly_set_fmpz_mod_bpoly(pfac->coeffs + i, A->bits, F->coeffs + i, 0, 1, ctx);
        fmpz_mod_mpoly_make_monic(pfac->coeffs + i, pfac->coeffs + i, ctx);
    }

    /* number of of local factors can only decrease from here on */
    subset = FLINT_ARRAY_REALLOC(subset, pfac->length, slong);

    for (m = 2; m <= n; m++)
    {
        fmpz_mod_mpoly_set(q, m < n ? Aevals + m : A, ctx);
        fmpz_mod_mpoly_set(p, Aevals + m - 1, ctx);

    #if FLINT_WANT_ASSERT
        fmpz_mod_mpoly_one(t, ctx);
        for (i = 0; i < pfac->length; i++)
            fmpz_mod_mpoly_mul(t, t, pfac->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mod_mpoly_equal(t, p, ctx));
    #endif

        /* if one local factor, A must be irreducible */
        if (pfac->length < 2)
        {
            fmpz_mod_mpolyv_fit_length(fac, 1, ctx);
            fac->length = 1;
            fmpz_mod_mpoly_set(fac->coeffs + 0, A, ctx);
            success = 1;
            goto cleanup;
        }

        success = _try_lift(qfac, q, pfac, p, m, alphas, n, ctx);
        if (success > 0)
        {
            fmpz_mod_mpolyv_swap(qfac, pfac, ctx);
            continue;
        }
        else if (success < 0)
        {
            success = 0;
            goto cleanup;
        }

        /* if we couldn't lift two local, A must be irreducible */
        if (pfac->length == 2)
        {
            fmpz_mod_mpolyv_fit_length(fac, 1, ctx);
            fac->length = 1;
            fmpz_mod_mpoly_set(fac->coeffs + 0, A, ctx);
            success = 1;
            goto cleanup;
        }

        qfac->length = 0;

        len = pfac->length;
        for (k = 0; k < len; k++)
            subset[k] = k;

        for (k = 1; k <= len/2; k++)
        {
            zassenhaus_subset_first(subset, len, k);

        #if FLINT_WANT_ASSERT
            fmpz_mod_mpoly_one(t, ctx);
            for (i = 0; i < len; i++)
            {
                int in = subset[i] >= 0;
                fmpz_mod_mpoly_mul(t, t, pfac->coeffs +
                                       (in ? subset[i] : -1 - subset[i]), ctx);
            }
            FLINT_ASSERT(fmpz_mod_mpoly_equal(t, p, ctx));
        #endif

            while (1)
            {
                fmpz_mod_mpolyv_fit_length(dfac, 2, ctx);
                dfac->length = 2;
                fmpz_mod_mpoly_one(dfac->coeffs + 0, ctx);
                fmpz_mod_mpoly_one(dfac->coeffs + 1, ctx);
                for (i = 0; i < len; i++)
                {
                    int in = subset[i] >= 0;
                    fmpz_mod_mpoly_mul(dfac->coeffs + in, dfac->coeffs + in,
                        pfac->coeffs + (in ? subset[i] : -1 - subset[i]), ctx);
                }

                success = _try_lift(tfac, q, dfac, p, m, alphas, n, ctx);
                if (success > 0)
                {
                    fmpz_mod_mpolyv_fit_length(qfac, qfac->length + 1, ctx);
                    fmpz_mod_mpoly_swap(qfac->coeffs + qfac->length,
                                                        tfac->coeffs + 1, ctx);
                    qfac->length++;

                    fmpz_mod_mpoly_swap(q, tfac->coeffs + 0, ctx);
                    fmpz_mod_mpoly_swap(p, dfac->coeffs + 0, ctx);
                    len -= k;
                    if (!zassenhaus_subset_next_disjoint(subset, len + k))
                        break;
                }
                else if (success < 0)
                {
                    success = 0;
                    goto cleanup;
                }
                else
                {
                    if (!zassenhaus_subset_next(subset, len))
                        break;
                }
            }
        }

        /* remnants are irreducible */
        if (!fmpz_mod_mpoly_is_fmpz(q, ctx))
        {
            fmpz_mod_mpolyv_fit_length(qfac, qfac->length + 1, ctx);
            fmpz_mod_mpoly_swap(qfac->coeffs + qfac->length, q, ctx);
            qfac->length++;
        }
        else
        {
            FLINT_ASSERT(fmpz_mod_mpoly_is_one(q, ctx));
        }

        fmpz_mod_mpolyv_swap(qfac, pfac, ctx);
    }

    success = 1;

    fmpz_mod_mpolyv_swap(fac, pfac, ctx);

cleanup: 

    flint_free(subset);
    _fmpz_vec_clear(alphas, n);
    for (i = 0; i < n; i++)
        fmpz_mod_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);
    flint_free(deg);
    flint_free(degeval);
    fmpz_mod_mpolyv_clear(pfac, ctx);
    fmpz_mod_mpolyv_clear(qfac, ctx);
    fmpz_mod_mpolyv_clear(tfac, ctx);
    fmpz_mod_mpolyv_clear(dfac, ctx);
    fmpz_mod_mpoly_clear(t, ctx);
    fmpz_mod_mpoly_clear(p, ctx);
    fmpz_mod_mpoly_clear(q, ctx);
    fmpz_mod_poly_clear(c, ctx->ffinfo);
    fmpz_mod_bpoly_clear(B, ctx->ffinfo);
    fmpz_mod_tpoly_clear(F, ctx->ffinfo);

#if FLINT_WANT_ASSERT
    if (success)
    {
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_one(t, ctx);
        for (i = 0; i < fac->length; i++)
            fmpz_mod_mpoly_mul(t, t, fac->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mod_mpoly_equal(t, A, ctx));
        fmpz_mod_mpoly_clear(t, ctx);
    }
#endif

    return success;
}
