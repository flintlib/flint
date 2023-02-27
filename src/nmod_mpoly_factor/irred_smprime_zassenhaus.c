/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "fq_nmod_mpoly_factor.h"

/*
    return:
        1: success
        0: lift is impossible
       -1: failed, don't try again
*/
static int _try_lift(
    nmod_mpolyv_t qfac,
    const nmod_mpoly_t q,
    const nmod_mpolyv_t pfac,
    const nmod_mpoly_t p,
    slong m,
    mp_limb_t * alpha,
    slong n,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    slong * newdeg;
    nmod_mpoly_t lcq, lcp, t, newq;

    FLINT_ASSERT(pfac->length > 1);

    newdeg = (slong *) flint_malloc((n + 1)*sizeof(slong));
    nmod_mpoly_init(lcq, ctx);
    nmod_mpoly_init(lcp, ctx);
    nmod_mpoly_init(t, ctx);
    nmod_mpoly_init(newq, ctx);

#if FLINT_WANT_ASSERT
    nmod_mpoly_one(t, ctx);
    for (i = 0; i < pfac->length; i++)
        nmod_mpoly_mul(t, t, pfac->coeffs + i, ctx);
    FLINT_ASSERT(nmod_mpoly_equal(t, p, ctx));
#endif

    _nmod_mpoly_get_lead0(lcq, q, ctx);
    nmod_mpoly_evaluate_one_ui(lcp, lcq, m, alpha[m - 1], ctx);

    FLINT_ASSERT(lcp->length > 0);

    nmod_mpoly_pow_ui(t, lcq, pfac->length - 1, ctx);
    nmod_mpoly_mul(newq, q, t, ctx);

    if (newq->bits > FLINT_BITS)
    {
        success = -1;
        goto cleanup;
    }

    nmod_mpoly_degrees_si(newdeg, newq, ctx);

    nmod_mpolyv_fit_length(qfac, pfac->length, ctx);
    qfac->length = pfac->length;
    for (i = 0; i < pfac->length; i++)
    {
        _nmod_mpoly_get_lead0(t, pfac->coeffs + i, ctx);
        success = nmod_mpoly_divides(t, lcp, t, ctx);
        FLINT_ASSERT(success);
        nmod_mpoly_mul(qfac->coeffs + i, pfac->coeffs + i, t, ctx);
        _nmod_mpoly_set_lead0(qfac->coeffs + i, qfac->coeffs + i, lcq, ctx);
    }

    success = nmod_mpoly_hlift(m, qfac->coeffs, qfac->length,
                                                     alpha, newq, newdeg, ctx);
    if (!success)
        goto cleanup;

    for (i = 0; i < qfac->length; i++)
    {
        /* hlift should not have returned any large bits */
        FLINT_ASSERT(qfac->coeffs[i].bits <= FLINT_BITS);

        if (!nmod_mpolyl_content(t, qfac->coeffs + i, 1, ctx))
        {
            success = -1;
            goto cleanup;
        }

        success = nmod_mpoly_divides(qfac->coeffs + i, qfac->coeffs + i, t, ctx);
        FLINT_ASSERT(success);
        nmod_mpoly_make_monic(qfac->coeffs + i, qfac->coeffs + i, ctx);
    }

    success = 1;

cleanup:

    flint_free(newdeg);
    nmod_mpoly_clear(lcq, ctx);
    nmod_mpoly_clear(lcp, ctx);
    nmod_mpoly_clear(t, ctx);
    nmod_mpoly_clear(newq, ctx);

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        nmod_mpoly_init(t, ctx);
        nmod_mpoly_one(t, ctx);
        for (i = 0; i < qfac->length; i++)
            nmod_mpoly_mul(t, t, qfac->coeffs + i, ctx);
        FLINT_ASSERT(nmod_mpoly_equal(t, q, ctx));
        nmod_mpoly_clear(t, ctx);
    }
#endif

    return success;
}


int nmod_mpoly_factor_irred_smprime_zassenhaus(
    nmod_mpolyv_t fac,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    int tries_remaining = 10;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, m, len;
    slong * subset;
    mp_limb_t * alpha;
    nmod_mpoly_struct * Aevals;
    slong * deg, * degeval;
    nmod_mpolyv_t qfac, pfac, tfac, dfac;
    nmod_mpoly_t t, p, q;
    n_poly_t c;
    n_bpoly_t B;
    n_tpoly_t F;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 1);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    subset = (slong*) flint_malloc(4*sizeof(slong));
    alpha = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));
    Aevals    = (nmod_mpoly_struct *) flint_malloc(n*sizeof(nmod_mpoly_struct));
    deg     = (slong *) flint_malloc((n + 1)*sizeof(slong));
    degeval = (slong *) flint_malloc((n + 1)*sizeof(slong));
    for (i = 0; i < n; i++)
        nmod_mpoly_init(Aevals + i, ctx);
    nmod_mpolyv_init(pfac, ctx);
    nmod_mpolyv_init(qfac, ctx);
    nmod_mpolyv_init(tfac, ctx);
    nmod_mpolyv_init(dfac, ctx);
    nmod_mpoly_init(t, ctx);
    nmod_mpoly_init(p, ctx);
    nmod_mpoly_init(q, ctx);
    n_poly_init(c);
    n_bpoly_init(B);
    n_tpoly_init(F);

    nmod_mpoly_degrees_si(deg, A, ctx);

next_alpha:

    if (--tries_remaining < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < n; i++)
        alpha[i] = n_urandint(state, ctx->mod.n - 1) + 1;

    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        nmod_mpoly_evaluate_one_ui(Aevals + i, i == n - 1 ? A :
                                         Aevals + i + 1, i + 1, alpha[i], ctx);
        nmod_mpoly_repack_bits_inplace(Aevals + i, A->bits, ctx);
        nmod_mpoly_degrees_si(degeval, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
            if (degeval[j] != deg[j])
                goto next_alpha;
    }

    /* make sure univar is squarefree */
    nmod_mpoly_derivative(t, Aevals + 0, 0, ctx);
    success = nmod_mpoly_gcd(t, t, Aevals + 0, ctx);
    if (!success)
        goto cleanup;
    if (!nmod_mpoly_is_one(t, ctx))
        goto next_alpha;

    /* make evaluations primitive */
    for (i = n - 1; i > 0; i--)
    {
        success = nmod_mpolyl_content(t, Aevals + i, 1, ctx);
        if (!success)
            goto cleanup;
        success = nmod_mpoly_divides(Aevals + i, Aevals + i, t, ctx);
        FLINT_ASSERT(success);
        nmod_mpoly_make_monic(Aevals + i, Aevals + i, ctx);
        nmod_mpoly_repack_bits_inplace(Aevals + i, A->bits, ctx);
    }

    nmod_mpoly_get_bpoly(B, Aevals + 1, 0, 1, ctx);
    success = n_bpoly_mod_factor_smprime(c, F, B, 1, ctx->mod);
    if (!success)
        goto next_alpha;

    FLINT_ASSERT(n_poly_degree(c) == 0);

    nmod_mpolyv_fit_length(pfac, F->length, ctx);
    pfac->length = F->length;
    for (i = 0; i < F->length; i++)
    {
        nmod_mpoly_set_bpoly(pfac->coeffs + i, A->bits, F->coeffs + i, 0, 1, ctx);
        nmod_mpoly_make_monic(pfac->coeffs + i, pfac->coeffs + i, ctx);
    }

    /* number of of local factors can only decrease from here on */
    subset = flint_realloc(subset, pfac->length*sizeof(slong));

    for (m = 2; m <= n; m++)
    {
        nmod_mpoly_set(q, m < n ? Aevals + m : A, ctx);
        nmod_mpoly_set(p, Aevals + m - 1, ctx);

    #if FLINT_WANT_ASSERT
        nmod_mpoly_one(t, ctx);
        for (i = 0; i < pfac->length; i++)
            nmod_mpoly_mul(t, t, pfac->coeffs + i, ctx);
        FLINT_ASSERT(nmod_mpoly_equal(t, p, ctx));
    #endif

        /* if one local factor, A must be irreducible */
        if (pfac->length < 2)
        {
            nmod_mpolyv_fit_length(fac, 1, ctx);
            fac->length = 1;
            nmod_mpoly_set(fac->coeffs + 0, A, ctx);
            success = 1;
            goto cleanup;
        }

        success = _try_lift(qfac, q, pfac, p, m, alpha, n, ctx);
        if (success > 0)
        {
            nmod_mpolyv_swap(qfac, pfac, ctx);
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
            nmod_mpolyv_fit_length(fac, 1, ctx);
            fac->length = 1;
            nmod_mpoly_set(fac->coeffs + 0, A, ctx);
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
            nmod_mpoly_one(t, ctx);
            for (i = 0; i < len; i++)
            {
                int in = subset[i] >= 0;
                nmod_mpoly_mul(t, t, pfac->coeffs +
                                       (in ? subset[i] : -1 - subset[i]), ctx);
            }
            FLINT_ASSERT(nmod_mpoly_equal(t, p, ctx));
        #endif

            while (1)
            {
                nmod_mpolyv_fit_length(dfac, 2, ctx);
                dfac->length = 2;
                nmod_mpoly_one(dfac->coeffs + 0, ctx);
                nmod_mpoly_one(dfac->coeffs + 1, ctx);
                for (i = 0; i < len; i++)
                {
                    int in = subset[i] >= 0;
                    nmod_mpoly_mul(dfac->coeffs + in, dfac->coeffs + in,
                        pfac->coeffs + (in ? subset[i] : -1 - subset[i]), ctx);
                }

                success = _try_lift(tfac, q, dfac, p, m, alpha, n, ctx);
                if (success > 0)
                {
                    nmod_mpolyv_fit_length(qfac, qfac->length + 1, ctx);
                    nmod_mpoly_swap(qfac->coeffs + qfac->length,
                                                        tfac->coeffs + 1, ctx);
                    qfac->length++;

                    nmod_mpoly_swap(q, tfac->coeffs + 0, ctx);
                    nmod_mpoly_swap(p, dfac->coeffs + 0, ctx);
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
        if (!nmod_mpoly_is_ui(q, ctx))
        {
            nmod_mpolyv_fit_length(qfac, qfac->length + 1, ctx);
            nmod_mpoly_swap(qfac->coeffs + qfac->length, q, ctx);
            qfac->length++;
        }
        else
        {
            FLINT_ASSERT(nmod_mpoly_is_one(q, ctx));
        }

        nmod_mpolyv_swap(qfac, pfac, ctx);
    }

    success = 1;

    nmod_mpolyv_swap(fac, pfac, ctx);

cleanup: 

    flint_free(subset);
    flint_free(alpha);
    for (i = 0; i < n; i++)
        nmod_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);
    flint_free(deg);
    flint_free(degeval);
    nmod_mpolyv_clear(pfac, ctx);
    nmod_mpolyv_clear(qfac, ctx);
    nmod_mpolyv_clear(tfac, ctx);
    nmod_mpolyv_clear(dfac, ctx);
    nmod_mpoly_clear(t, ctx);
    nmod_mpoly_clear(p, ctx);
    nmod_mpoly_clear(q, ctx);
    n_poly_clear(c);
    n_bpoly_clear(B);
    n_tpoly_clear(F);

#if FLINT_WANT_ASSERT
    if (success)
    {
        nmod_mpoly_init(t, ctx);
        nmod_mpoly_one(t, ctx);
        for (i = 0; i < fac->length; i++)
            nmod_mpoly_mul(t, t, fac->coeffs + i, ctx);
        FLINT_ASSERT(nmod_mpoly_equal(t, A, ctx));
        nmod_mpoly_clear(t, ctx);
    }
#endif

    return success;
}
