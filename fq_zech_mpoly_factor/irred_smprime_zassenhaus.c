/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"

/*
    return:
        1: success
        0: lift is impossible
       -1: failed, don't try again
*/
static int _try_lift(
    fq_zech_mpolyv_t qfac,
    const fq_zech_mpoly_t q,
    const fq_zech_mpolyv_t pfac,
    const fq_zech_mpoly_t p,
    slong m,
    fq_zech_struct * alpha,
    slong n,
    const fq_zech_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    slong * newdeg;
    fq_zech_mpoly_t lcq, lcp, t, newq;
    fq_zech_mpoly_univar_t u;

    FLINT_ASSERT(pfac->length > 1);

    newdeg = (slong *) flint_malloc((n + 1)*sizeof(slong));
    fq_zech_mpoly_init(lcq, ctx);
    fq_zech_mpoly_init(lcp, ctx);
    fq_zech_mpoly_init(t, ctx);
    fq_zech_mpoly_init(newq, ctx);
    fq_zech_mpoly_univar_init(u, ctx);

#if FLINT_WANT_ASSERT
    fq_zech_mpoly_one(t, ctx);
    for (i = 0; i < pfac->length; i++)
        fq_zech_mpoly_mul(t, t, pfac->coeffs + i, ctx);
    FLINT_ASSERT(fq_zech_mpoly_equal(t, p, ctx));
#endif

    _fq_zech_mpoly_get_lead0(lcq, q, ctx);
    fq_zech_mpoly_evaluate_one_fq_zech(lcp, lcq, m, alpha + m - 1, ctx);

    FLINT_ASSERT(lcp->length > 0);

    fq_zech_mpoly_pow_ui(t, lcq, pfac->length - 1, ctx);
    fq_zech_mpoly_mul(newq, q, t, ctx);

    if (newq->bits > FLINT_BITS)
    {
        success = -1;
        goto cleanup;
    }

    fq_zech_mpoly_degrees_si(newdeg, newq, ctx);

    fq_zech_mpolyv_fit_length(qfac, pfac->length, ctx);
    qfac->length = pfac->length;
    for (i = 0; i < pfac->length; i++)
    {
        _fq_zech_mpoly_get_lead0(t, pfac->coeffs + i, ctx);
        success = fq_zech_mpoly_divides(t, lcp, t, ctx);
        FLINT_ASSERT(success);
        fq_zech_mpoly_mul(qfac->coeffs + i, pfac->coeffs + i, t, ctx);
        _fq_zech_mpoly_set_lead0(qfac->coeffs + i, qfac->coeffs + i, lcq, ctx);
    }

    success = fq_zech_mpoly_hlift(m, qfac->coeffs, qfac->length,
                                                     alpha, newq, newdeg, ctx);
    if (!success)
        goto cleanup;

    for (i = 0; i < qfac->length; i++)
    {
        fq_zech_mpoly_to_univar(u, qfac->coeffs + i, 0, ctx);
        success = fq_zech_mpoly_univar_content_mpoly(t, u, ctx);
        if (!success)
        {
            success = -1;
            goto cleanup;
        }
        success = fq_zech_mpoly_divides(qfac->coeffs + i,
                                        qfac->coeffs + i, t, ctx);
        FLINT_ASSERT(success);
        fq_zech_mpoly_make_monic(qfac->coeffs + i, qfac->coeffs + i, ctx);
    }

    success = 1;

cleanup:

    flint_free(newdeg);
    fq_zech_mpoly_clear(lcq, ctx);
    fq_zech_mpoly_clear(lcp, ctx);
    fq_zech_mpoly_clear(t, ctx);
    fq_zech_mpoly_clear(newq, ctx);
    fq_zech_mpoly_univar_clear(u, ctx);

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        fq_zech_mpoly_init(t, ctx);
        fq_zech_mpoly_one(t, ctx);
        for (i = 0; i < qfac->length; i++)
            fq_zech_mpoly_mul(t, t, qfac->coeffs + i, ctx);
        FLINT_ASSERT(fq_zech_mpoly_equal(t, q, ctx));
        fq_zech_mpoly_clear(t, ctx);
    }
#endif

    return success;
}

/*
    return:
        1: success
        0: failed, may try again
       -1: failed, don't try again
*/
int fq_zech_mpoly_factor_irred_smprime_zassenhaus(
    fq_zech_mpolyv_t fac,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int tries_left = 10;
    int success;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, m, len;
    slong * subset;
    fq_zech_struct * alpha;
    fq_zech_mpoly_struct * Aevals;
    slong * deg, * degeval;
    fq_zech_mpolyv_t qfac, pfac, tfac, dfac;
    fq_zech_mpoly_t t, p, q;
    fq_zech_mpoly_univar_t u;
    fq_zech_poly_t c;
    fq_zech_bpoly_t B;
    fq_zech_tpoly_t F;

    subset = (slong*) flint_malloc(4*sizeof(slong));
	alpha = FLINT_ARRAY_ALLOC(n, fq_zech_struct);
    for (i = 0; i < n; i++)
        fq_zech_init(alpha + i, ctx->fqctx);

    deg     = FLINT_ARRAY_ALLOC(n + 1, slong);
    degeval = FLINT_ARRAY_ALLOC(n + 1, slong);

    Aevals = FLINT_ARRAY_ALLOC(n, fq_zech_mpoly_struct);
	for (i = 0; i < n; i++)
		fq_zech_mpoly_init(Aevals + i, ctx);

    fq_zech_mpolyv_init(pfac, ctx);
	fq_zech_mpolyv_init(qfac, ctx);
    fq_zech_mpolyv_init(tfac, ctx);
	fq_zech_mpolyv_init(dfac, ctx);
	fq_zech_mpoly_init(t, ctx);
	fq_zech_mpoly_init(p, ctx);
	fq_zech_mpoly_init(q, ctx);
	fq_zech_mpoly_univar_init(u, ctx);
    fq_zech_poly_init(c, ctx->fqctx);
    fq_zech_bpoly_init(B, ctx->fqctx);
    fq_zech_tpoly_init(F, ctx->fqctx);

	fq_zech_mpoly_degrees_si(deg, A, ctx);
    goto got_alpha;

next_alpha:

    tries_left--;
    if (tries_left < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < n; i++)
        fq_zech_rand(alpha + i, state, ctx->fqctx);

got_alpha:

	/* ensure degrees do not drop under evalutaion */
    for (i = n - 1; i >= 0; i--)
    {
    	fq_zech_mpoly_evaluate_one_fq_zech(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
    	fq_zech_mpoly_degrees_si(degeval, Aevals + i, ctx);
	    for (j = 0; j <= i; j++)
        {
	    	if (degeval[j] != deg[j])
			    goto next_alpha;
        }
    }

	fq_zech_mpoly_get_fq_zech_poly(c, Aevals + 0, 0, ctx);
	if (!fq_zech_poly_is_squarefree(c, ctx->fqctx))
		goto next_alpha;

    /* make evaluations primitive */
    for (i = n - 1; i > 0; i--)
    {
        fq_zech_mpoly_to_univar(u, Aevals + i, 0, ctx);
        success = fq_zech_mpoly_univar_content_mpoly(t, u, ctx);
        if (!success)
            goto cleanup;
        success = fq_zech_mpoly_divides(Aevals + i, Aevals + i, t, ctx);
        FLINT_ASSERT(success);
        fq_zech_mpoly_make_monic(Aevals + i, Aevals + i, ctx);
    }

    fq_zech_mpoly_get_fq_zech_bpoly(B, Aevals + 1, 0, 1, ctx);
    success = fq_zech_bpoly_factor_smprime(c, F, B, 1, ctx->fqctx);
    if (!success)
        goto next_alpha;

    FLINT_ASSERT(fq_zech_poly_degree(c, ctx->fqctx) == 0);

    fq_zech_mpolyv_fit_length(pfac, F->length, ctx);
    pfac->length = F->length;
    for (i = 0; i < F->length; i++)
    {
        fq_zech_mpoly_set_fq_zech_bpoly(pfac->coeffs + i, A->bits,
                                                     F->coeffs + i, 0, 1, ctx);
        fq_zech_mpoly_make_monic(pfac->coeffs + i, pfac->coeffs + i, ctx);
    }

    for (m = 2; m <= n; m++)
    {
        fq_zech_mpoly_set(q, m < n ? Aevals + m : A, ctx);
        fq_zech_mpoly_set(p, Aevals + m - 1, ctx);

    #if FLINT_WANT_ASSERT
        fq_zech_mpoly_one(t, ctx);
        for (i = 0; i < pfac->length; i++)
            fq_zech_mpoly_mul(t, t, pfac->coeffs + i, ctx);
        FLINT_ASSERT(fq_zech_mpoly_equal(t, p, ctx));
    #endif

        if (pfac->length < 2)
        {
	        fq_zech_mpolyv_fit_length(fac, 1, ctx);
            fac->length = 1;
            fq_zech_mpoly_set(fac->coeffs + 0, A, ctx);
		    success = 1;
		    goto cleanup;
        }

        success = _try_lift(qfac, q, pfac, p, m, alpha, n, ctx);
        if (success > 0)
        {
            fq_zech_mpolyv_swap(qfac, pfac, ctx);
            continue;
        }
        else if (success < 0)
        {
            goto cleanup;
        }

        if (pfac->length == 2)
        {
	        fq_zech_mpolyv_fit_length(fac, 1, ctx);
            fac->length = 1;
            fq_zech_mpoly_set(fac->coeffs + 0, A, ctx);
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
            fq_zech_mpoly_one(t, ctx);
            for (i = 0; i < len; i++)
            {
                int in = subset[i] >= 0;
                fq_zech_mpoly_mul(t, t, pfac->coeffs +
                                       (in ? subset[i] : -1 - subset[i]), ctx);
            }
            FLINT_ASSERT(fq_zech_mpoly_equal(t, p, ctx));
        #endif

            while (1)
            {
                fq_zech_mpolyv_fit_length(dfac, 2, ctx);
                dfac->length = 2;
                fq_zech_mpoly_one(dfac->coeffs + 0, ctx);
                fq_zech_mpoly_one(dfac->coeffs + 1, ctx);
                for (i = 0; i < len; i++)
                {
                    int in = subset[i] >= 0;
                    fq_zech_mpoly_mul(dfac->coeffs + in, dfac->coeffs + in,
                        pfac->coeffs + (in ? subset[i] : -1 - subset[i]), ctx);
                }

                success = _try_lift(tfac, q, dfac, p, m, alpha, n, ctx);
                if (success > 0)
                {
                    fq_zech_mpolyv_fit_length(qfac, qfac->length + 1, ctx);
                    fq_zech_mpoly_swap(qfac->coeffs + qfac->length,
                                                        tfac->coeffs + 1, ctx);
                    qfac->length++;

                    fq_zech_mpoly_swap(q, tfac->coeffs + 0, ctx);
                    fq_zech_mpoly_swap(p, dfac->coeffs + 0, ctx);
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
        if (!fq_zech_mpoly_is_fq_zech(q, ctx))
        {
            fq_zech_mpolyv_fit_length(qfac, qfac->length + 1, ctx);
            fq_zech_mpoly_swap(qfac->coeffs + qfac->length, q, ctx);
            qfac->length++;
        }
        else
        {
            FLINT_ASSERT(fq_zech_mpoly_is_one(q, ctx));
        }

        fq_zech_mpolyv_swap(qfac, pfac, ctx);  
    }

    success = 1;

    fq_zech_mpolyv_swap(fac, pfac, ctx);

cleanup:

    flint_free(subset);

    for (i = 0; i < n; i++)
        fq_zech_clear(alpha + i, ctx->fqctx);
    flint_free(alpha);

    for (i = 0; i < n; i++)
        fq_zech_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);

    flint_free(deg);
    flint_free(degeval);

    fq_zech_mpolyv_clear(pfac, ctx);
    fq_zech_mpolyv_clear(qfac, ctx);
    fq_zech_mpolyv_clear(tfac, ctx);
    fq_zech_mpolyv_clear(dfac, ctx);
    fq_zech_mpoly_clear(t, ctx);
    fq_zech_mpoly_clear(p, ctx);
    fq_zech_mpoly_clear(q, ctx);
    fq_zech_mpoly_univar_clear(u, ctx);
    fq_zech_poly_clear(c, ctx->fqctx);
    fq_zech_bpoly_clear(B, ctx->fqctx);
    fq_zech_tpoly_clear(F, ctx->fqctx);

    FLINT_ASSERT(success == 0 || success == 1);

#if FLINT_WANT_ASSERT
    if (success)
    {
        fq_zech_mpoly_init(t, ctx);
        fq_zech_mpoly_one(t, ctx);
        for (i = 0; i < fac->length; i++)
            fq_zech_mpoly_mul(t, t, fac->coeffs + i, ctx);
        FLINT_ASSERT(fq_zech_mpoly_equal(t, A, ctx));
        fq_zech_mpoly_clear(t, ctx);
    }
#endif

    return success;
}

