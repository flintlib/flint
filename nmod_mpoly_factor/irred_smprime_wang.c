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
    only E and alphas are shifted by "var"
    so output is in E[0]
*/
int _nmod_mpoly_evaluate_rest_n_poly(
    n_poly_struct * E,
    slong * starts,
    slong * ends,
    slong * stops,
    ulong * es,
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    slong var,
    const n_poly_struct * alphas,
    const slong * offsets,
    const slong * shifts,
    slong N,
    ulong mask,
    slong nvars,
    nmod_t ctx)
{
    slong v, stop;
    ulong next_e;

    FLINT_ASSERT(var < nvars);

    E -= var;
    alphas -= var;

    v = var;
    starts[v] = 0;
    ends[v] = Alen;
    n_poly_zero(E + v);

    if (Alen < 1)
        return 1;

calculate:
/*
    input:
        v
        starts[v]
        ends[v]
*/
    FLINT_ASSERT(ends[v] > starts[v]);
    es[v] = mask & (Aexps[N*starts[v] + offsets[v]] >> shifts[v]);

    n_poly_zero(E + v);

next:

    FLINT_ASSERT(es[v] == (mask & (Aexps[N*starts[v] + offsets[v]] >> shifts[v])));

    stop = starts[v] + 1;
    while (stop < ends[v] &&
           (mask & (Aexps[N*stop + offsets[v]] >> shifts[v])) == es[v])
    {
        stop++;
    }
    stops[v] = stop;

    if (v + 1 < nvars)
    {
        starts[v + 1] = starts[v];
        ends[v + 1] = stops[v];
        v++;
        goto calculate;
calculate_return:
        n_poly_mod_add(E + v, E + v, E + v + 1, ctx);
    }
    else
    {
        n_poly_mod_add_ui(E + v, E + v, Acoeffs[starts[v]], ctx);
    }

    if (stops[v] < ends[v])
    {
        next_e = mask & (Aexps[N*stops[v] + offsets[v]] >> shifts[v]);
        FLINT_ASSERT(next_e < es[v]);
        n_poly_mod_pow(E + v + 1, alphas + v, es[v] - next_e, ctx);
        n_poly_mod_mul(E + v, E + v, E + v + 1, ctx);
        es[v] = next_e;
        starts[v] = stops[v];
        goto next;
    }
    else
    {
        n_poly_mod_pow(E + v + 1, alphas + v, es[v], ctx);
        n_poly_mod_mul(E + v, E + v, E + v + 1, ctx);
    }

    if (v > var)
    {
        v--;
        goto calculate_return;
    }

    return 1;
}


void _eval_to_bpoly(
    n_bpoly_t E,
    const nmod_mpoly_t A,
    const n_poly_struct * alphabetas,
    const nmod_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong * offsets, * shifts;
    slong offset, shift;
    slong start, stop;
    ulong e, mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong * starts, * ends, * stops;
    ulong * es;
    n_poly_struct * realE;

    E->length = 0;
    if (A->length < 1)
        return;

    starts = FLINT_ARRAY_ALLOC(n, slong);
    ends   = FLINT_ARRAY_ALLOC(n, slong);
    stops  = FLINT_ARRAY_ALLOC(n, slong);
    es     = FLINT_ARRAY_ALLOC(n, ulong);
    realE  = FLINT_ARRAY_ALLOC(n + 1, n_poly_struct);
    for (i = 0; i < n + 1; i++)
        n_poly_init(realE + i);

    offsets = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    for (i = 0; i < ctx->minfo->nvars; i++)
        mpoly_gen_offset_shift_sp(offsets + i, shifts + i, i, A->bits, ctx->minfo);

    offset = offsets[0];
    shift = shifts[0];

    start = 0;
    e = mask & (A->exps[N*start + offset] >> shift);

next:

    FLINT_ASSERT(start < A->length);
    FLINT_ASSERT(e == (mask & (A->exps[N*start + offset] >> shift)));

    stop = start + 1;
    while (stop < A->length && (mask & (A->exps[N*stop + offset] >> shift)) == e)
        stop++;

    n_bpoly_fit_length(E, e + 1);
    while (E->length <= e)
    {
        n_poly_zero(E->coeffs + E->length);
        E->length++;
    }

    _nmod_mpoly_evaluate_rest_n_poly(realE, starts, ends, stops, es,
                    A->coeffs + start, A->exps + N*start, stop - start, 1,
                                        alphabetas, offsets, shifts, N, mask,
                                          ctx->minfo->nvars, ctx->ffinfo->mod);
    n_poly_set(E->coeffs + e, realE + 0);

    if (stop < A->length)
    {
        FLINT_ASSERT(e > (mask & (A->exps[N*stop + offset] >> shift)));
        e = (mask & (A->exps[N*stop + offset] >> shift));
        start = stop;
        goto next;
    }

    n_bpoly_normalise(E);

    for (i = 0; i < n + 1; i++)
        n_poly_clear(realE + i);
    flint_free(realE);
    flint_free(starts);
    flint_free(ends);
    flint_free(stops);
    flint_free(es);

    flint_free(offsets);
    flint_free(shifts);
}


/* A = B(gen(var), 0) */
void _nmod_mpoly_set_bpoly_var1_zero(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const n_bpoly_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);
    slong i, Alen;
    slong Blen = B->length;
    ulong * genexp;
    TMP_INIT;

    TMP_START;

    genexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    if (Abits <= FLINT_BITS)
        mpoly_gen_monomial_sp(genexp, var, Abits, ctx->minfo);
    else
        mpoly_gen_monomial_offset_mp(genexp, var, Abits, ctx->minfo);

    Alen = 2;
    for (i = 0; i < Blen; i++)
        Alen += (B->coeffs[i].length > 0);

    nmod_mpoly_fit_length_set_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        mp_limb_t c = n_poly_get_coeff(B->coeffs + i, 0);
        if (c == 0)
            continue;

        FLINT_ASSERT(Alen < A->alloc);
        A->coeffs[Alen] = c;
        if (Abits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps + N*Alen, genexp, N, i);
        else
            mpoly_monomial_mul_ui_mp(A->exps + N*Alen, genexp, N, i);
        Alen++;
    }
    A->length = Alen;

    TMP_END;
}

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
    {
        alpha[i] = n_urandint(state, ctx->ffinfo->mod.n - 1) + 1;
    }

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
	nmod_mpoly_gcd(t, t, Aevals + 0, ctx);
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
            alphabetas[i].coeffs[j] = n_urandint(state, ctx->ffinfo->mod.n);
        alphabetas[i].length = alphabetas_length;
        _n_poly_normalise(alphabetas + i);
    }

    _eval_to_bpoly(Ab, A, alphabetas, ctx);
    success = n_bpoly_mod_factor_smprime(Abfc, Abfp, Ab, 0, ctx->ffinfo->mod);
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
        _nmod_mpoly_set_bpoly_var1_zero(fac->coeffs + i, newA->bits, Abfp->coeffs + i, 0, ctx);
        FLINT_ASSERT(fac->coeffs[i].length > 0);
        q = nmod_inv(fac->coeffs[i].coeffs[0], ctx->ffinfo->mod);
        q = nmod_mul(q, new_lcs->coeffs[0*r + i].coeffs[0], ctx->ffinfo->mod);
        nmod_mpoly_scalar_mul_nmod_invertible(fac->coeffs + i, fac->coeffs + i, q, ctx);
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
        nmod_mpoly_univar_t u;
        nmod_mpoly_univar_init(u, ctx);
        for (i = 0; i < r; i++)
        {
            nmod_mpoly_to_univar(u, fac->coeffs + i, 0, ctx);
            success = _nmod_mpoly_vec_content_mpoly(t, u->coeffs, u->length, ctx);
            if (!success)
            {
                nmod_mpoly_univar_clear(u, ctx);
                goto cleanup;
            }
            success = nmod_mpoly_divides(fac->coeffs + i,
                                         fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
        }
        nmod_mpoly_univar_clear(u, ctx);
    }

    for (i = 0; i < r; i++)
        nmod_mpoly_make_monic(fac->coeffs + i, fac->coeffs + i, ctx);

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

#if WANT_ASSERT
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
/*
flint_printf("nmod_mpoly_factor_irred_smprime_wang p = %wu returning %d\n", ctx->ffinfo->mod.n, success);
*/
	return success;
}
