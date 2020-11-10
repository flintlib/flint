/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"

#if FLINT_WANT_ASSERT
/*
    return: 0 no
            1 yes
           -1 don't know
*/
static int fq_nmod_mpoly_factor_is_pairwise_prime(
	const fq_nmod_mpoly_factor_t f, 
	const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
	int result;
	slong i, j;
	fq_nmod_mpoly_t g;

	fq_nmod_mpoly_init(g, ctx);

	for (i = 0; i + 1 < f->num; i++)
	for (j = i + 1; j < f->num; j++)
	{
		/* make sure factors are monic */
		if (f->poly[i].length < 1 ||
            f->poly[j].length < 1 ||
            !_n_fq_is_one(f->poly[i].coeffs + d*0, d) ||
            !_n_fq_is_one(f->poly[j].coeffs + d*0, d))
		{
			result = 0;
			goto cleanup;
		}

		if (fq_nmod_mpoly_gcd(g, f->poly + i, f->poly + j, ctx))
		{
			if (!fq_nmod_mpoly_is_one(g, ctx))
			{
				result = 0;
				goto cleanup;
			}
		}
		else
		{
			result = -1;
			goto cleanup;
		}
	}

	result = 1;

cleanup:

	fq_nmod_mpoly_clear(g, ctx);

	return result;
}
#endif

/*
	b and c are pairwise prime
	produce a=b*c pairwise prime and return 1
		else return 0 with a undefined
*/
static int fq_nmod_mpoly_factor_mul_pairwise_prime(
	fq_nmod_mpoly_factor_t a,
	fq_nmod_mpoly_factor_t b,   /* clobbered */
	fq_nmod_mpoly_factor_t c,   /* clobbered */
	const fq_nmod_mpoly_ctx_t ctx)
{
	int success;
	slong i, j;
	fq_nmod_mpoly_t T1, T2;
	fq_nmod_mpoly_t g;
	fmpz_t t;
#if FLINT_WANT_ASSERT
    int ae_ok = 1;
	fq_nmod_mpoly_t ae, be, ce;
#endif

	if (a == b || a == c)
	{
		fq_nmod_mpoly_factor_t ta;
		fq_nmod_mpoly_factor_init(ta, ctx);
		success = fq_nmod_mpoly_factor_mul_pairwise_prime(ta, b, c, ctx);
		fq_nmod_mpoly_factor_swap(a, ta, ctx);
		fq_nmod_mpoly_factor_clear(ta, ctx);
		return success;
	}

#if FLINT_WANT_ASSERT
	fq_nmod_mpoly_init(ae, ctx);
	fq_nmod_mpoly_init(be, ctx);
	fq_nmod_mpoly_init(ce, ctx);

	ae_ok = fq_nmod_mpoly_factor_expand(be, b, ctx) &&
	        fq_nmod_mpoly_factor_expand(ce, c, ctx);
    if (ae_ok)
    	fq_nmod_mpoly_mul(ae, be, ce, ctx);
#endif

	fmpz_init(t);
	fq_nmod_mpoly_init(T1, ctx);
	fq_nmod_mpoly_init(T2, ctx);

	FLINT_ASSERT(fq_nmod_mpoly_factor_is_pairwise_prime(b, ctx) != 0);
	FLINT_ASSERT(fq_nmod_mpoly_factor_is_pairwise_prime(c, ctx) != 0);

    fq_nmod_mpoly_init(g, ctx);

	fq_nmod_mul(a->constant, b->constant, c->constant, ctx->fqctx);
	a->num = 0;
	for (i = 0; i < b->num; i++)
	for (j = 0; j < c->num; j++)
	{
		if (!fq_nmod_mpoly_gcd_cofactors(g, b->poly + i, c->poly + j, 
                                            b->poly + i, c->poly + j, ctx))
		{
			success = 0;
			goto cleanup;
		}

        if (!fq_nmod_mpoly_is_one(g, ctx))
        {
            fq_nmod_mpoly_factor_fit_length(a, a->num + 1, ctx);
            fq_nmod_mpoly_swap(a->poly + a->num, g, ctx);
            fmpz_add(a->exp + a->num, b->exp + i, c->exp + j);
            a->num++;
        }
	}

	for (i = 0; i < b->num; i++)
	{
        if (!fq_nmod_mpoly_is_one(b->poly + i, ctx))
        {
            fq_nmod_mpoly_factor_fit_length(a, a->num + 1, ctx);
            fq_nmod_mpoly_swap(a->poly + a->num, b->poly + i, ctx);
            fmpz_swap(a->exp + a->num, b->exp + i);
            a->num++;
        }
	}

	for (j = 0; j < c->num; j++)
	{
        if (!fq_nmod_mpoly_is_one(c->poly + j, ctx))
        {
            fq_nmod_mpoly_factor_fit_length(a, a->num + 1, ctx);
            fq_nmod_mpoly_swap(a->poly + a->num, c->poly + j, ctx);
            fmpz_swap(a->exp + a->num, c->exp + j);
            a->num++;
        }
	}

	success = 1;

cleanup:

    FLINT_ASSERT(!success || fq_nmod_mpoly_factor_is_pairwise_prime(a, ctx) != 0);
    FLINT_ASSERT(!(ae_ok && success) || fq_nmod_mpoly_factor_matches(ae, a, ctx));

#if FLINT_WANT_ASSERT
	fq_nmod_mpoly_clear(ae, ctx);
	fq_nmod_mpoly_clear(be, ctx);
	fq_nmod_mpoly_clear(ce, ctx);
#endif

    fq_nmod_mpoly_clear(g, ctx);

	fq_nmod_mpoly_clear(T1, ctx);
	fq_nmod_mpoly_clear(T2, ctx);
	fmpz_clear(t);

	return success;
}


/*
    a has zero derivative wrt gen(i) for all i for which vars_left[i] = 0
    gcd(a, derivative(a, gen(var))) = 1

    sep = false: just tack on a
    sep = true: make sure either
                1. derivative(a, gen(var)) = 0, or
                2. gcd(a, derivative(a, gen(i))) = 1
            holds for the other factors for all i for which vars_left[i] != 0
*/
static int _append_factor_sep(
    fq_nmod_mpoly_factor_t f,
    fq_nmod_mpoly_t a,
    ulong k,
    int * vars_left,
    const fq_nmod_mpoly_ctx_t ctx,
    int sep,
    fq_nmod_mpoly_t t) /* temp */
{
    slong v, org = f->num;

    if (fq_nmod_mpoly_is_fq_nmod(a, ctx))
    {
        FLINT_ASSERT(fq_nmod_mpoly_is_one(a, ctx));
        return 1;
    }

    fq_nmod_mpoly_factor_fit_length(f, org + 1, ctx);
    fq_nmod_mpoly_swap(f->poly + org, a, ctx);
    fmpz_set_ui(f->exp + org, k);
    f->num = org + 1;

    if (!sep)
        return 1;

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        slong i = org;

        if (!vars_left[v])
            continue;

        while (i < f->num)
        {
            fq_nmod_mpoly_derivative(t, f->poly + i, v, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
            {
                /* f[i] has zero derivative */
                i++;
                continue;
            }

            fq_nmod_mpoly_factor_fit_length(f, f->num + 1, ctx);
            fmpz_set_ui(f->exp + f->num, k);

            if (!fq_nmod_mpoly_gcd_cofactors(f->poly + f->num, f->poly + i, t,
                                                          f->poly + i, t, ctx))
            {
                return 0;
            }

            if (fq_nmod_mpoly_is_fq_nmod(f->poly + f->num, ctx))
            {
                /* f[i] is comprime with its derivative */
                i++;
            }
            else
            {
                /* f[i] and f[end] at least got smaller */
                f->num++;
            }
        }
    }

    return 1;
}

int _fq_nmod_mpoly_factor_separable(
    fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    int sep)
{
	int success;
	slong i, j, v, var;
    ulong k;
    int * vars_left;
	fmpz * shift, * stride;
	fq_nmod_mpoly_factor_t Tf;
	fmpz_t g, gr, p, pk;
    fq_nmod_mpoly_t B, C, U, V, W, G;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(n_fq_is_one(A->coeffs + 0, ctx->fqctx));

    fq_nmod_mpoly_init(B, ctx);
    fq_nmod_mpoly_init(C, ctx);
    fq_nmod_mpoly_init(U, ctx);
    fq_nmod_mpoly_init(V, ctx);
    fq_nmod_mpoly_init(W, ctx);
    fq_nmod_mpoly_init(G, ctx);

	fmpz_init_set_ui(p, ctx->fqctx->modulus->mod.n);
	fmpz_init(pk);
	fmpz_init(g);
	fmpz_init(gr);
	fq_nmod_mpoly_factor_init(Tf, ctx);
	shift = _fmpz_vec_init(ctx->minfo->nvars);
	stride = _fmpz_vec_init(ctx->minfo->nvars);
    vars_left = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, int);
    for (v = 0; v < ctx->minfo->nvars; v++)
        vars_left[v] = 1;

    fq_nmod_mpoly_factor_one(f, ctx);
    fq_nmod_mpoly_set(C, A, ctx);

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        /* take the next variable with shortest derivative */
        var = -1;
        for (v = 0; v < ctx->minfo->nvars; v++)
        {
            if (!vars_left[v])
                continue;

            fq_nmod_mpoly_derivative(U, C, v, ctx);

            if (var < 0 || U->length < G->length)
            {
                var = v;
                fq_nmod_mpoly_swap(G, U, ctx);
            }
        }

        FLINT_ASSERT(var >= 0);
        FLINT_ASSERT(vars_left[var] == 1);

        vars_left[var] = 0;

        success = fq_nmod_mpoly_gcd_cofactors(C, W, V, C, G, ctx);
        if (!success)
            goto cleanup;

        for (k = 1; k + 1 < ctx->fqctx->modulus->mod.n &&
                            !(fq_nmod_mpoly_derivative(G, W, var, ctx),
                              fq_nmod_mpoly_sub(U, V, G, ctx),
                              fq_nmod_mpoly_is_zero(U, ctx)); k++)
        {
            success = fq_nmod_mpoly_gcd_cofactors(G, W, V, W, U, ctx);
            if (!success)
	            goto cleanup;

            success = _append_factor_sep(f, G, k, vars_left, ctx, sep, U);
            if (!success)
                goto cleanup;

            if (!fq_nmod_mpoly_is_one(W, ctx))
            {
                success = fq_nmod_mpoly_divides(U, C, W, ctx);
                FLINT_ASSERT(success);
                fq_nmod_mpoly_swap(C, U, ctx);
            }
        }

        success = _append_factor_sep(f, W, k, vars_left, ctx, sep, U);
        if (!success)
            goto cleanup;
	}

	if (fq_nmod_mpoly_is_fq_nmod(C, ctx))
	{
        FLINT_ASSERT(fq_nmod_mpoly_is_one(C, ctx));
	}
	else
	{
        slong d = fq_nmod_ctx_degree(ctx->fqctx);
        slong mk_mod_deg;

		fq_nmod_mpoly_deflation(shift, stride, C, ctx);
		fmpz_zero(g);
		for (var = 0; var < ctx->minfo->nvars; var++)
		{
			fmpz_gcd(g, g, stride + var);
			fmpz_gcd(g, g, shift + var);
		}

        k = fmpz_remove(gr, g, p);
        FLINT_ASSERT(k > 0);
		fmpz_pow_ui(pk, p, k);

        /* p^k th root of exponents */
		for (var = 0; var < ctx->minfo->nvars; var++)
		{
			fmpz_set(stride + var, pk);
			fmpz_zero(shift + var);
		}
		fq_nmod_mpoly_deflate(C, C, shift, stride, ctx);

        /* p^k th root of coefficients */
        mk_mod_deg = (ulong)k % (ulong)fq_nmod_ctx_degree(ctx->fqctx);
        if (mk_mod_deg > 0)
            mk_mod_deg = fq_nmod_ctx_degree(ctx->fqctx) - mk_mod_deg;

        for (i = 0; i < C->length; i++)
        {
            for (j = 0; j < mk_mod_deg; j++)
            {
                n_fq_pow_ui(C->coeffs + d*i, C->coeffs + d*i,
                                    fq_nmod_ctx_mod(ctx->fqctx).n, ctx->fqctx);
            }
        }

		success = _fq_nmod_mpoly_factor_separable(Tf, C, ctx, sep);
		if (!success)
			goto cleanup;

        /* p^k power of factors */
        /* f = f * Cf */
        FLINT_ASSERT(fq_nmod_is_one(Tf->constant, ctx->fqctx));
	    _fmpz_vec_scalar_mul_fmpz(Tf->exp, Tf->exp, Tf->num, pk);
		fq_nmod_mpoly_factor_mul_pairwise_prime(f, f, Tf, ctx);
	}

	success = 1;

cleanup:

    fq_nmod_mpoly_clear(C, ctx);
    fq_nmod_mpoly_clear(U, ctx);
    fq_nmod_mpoly_clear(V, ctx);
    fq_nmod_mpoly_clear(W, ctx);
    fq_nmod_mpoly_clear(G, ctx);

	fmpz_clear(p);
	fmpz_clear(pk);
	fmpz_clear(g);
	fmpz_clear(gr);
	fq_nmod_mpoly_factor_clear(Tf, ctx);
	_fmpz_vec_clear(shift, ctx->minfo->nvars);
	_fmpz_vec_clear(stride, ctx->minfo->nvars);
    flint_free(vars_left);

	return success;
}


/*
    if sep = true, each returned factor should satisfy:
        (1) monic
        (2) primitive wrt each variable
        (3) for all i, derivative(a, gen(i)) = 0, or
                       gcd(a, derivative(a, gen(i))) = 1
        (4) there is at least i for which derivative(a, gen(i)) != 0

    otherwise, the factors are just squarefree
*/
int fq_nmod_mpoly_factor_separable(
    fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    int sep)
{
    int success;
    slong i, j;
    fq_nmod_mpoly_factor_t g, t;

    if (!fq_nmod_mpoly_factor_content(f, A, ctx))
        return 0;

    fq_nmod_mpoly_factor_init(g, ctx);
    fq_nmod_mpoly_factor_init(t, ctx);

    fq_nmod_set(g->constant, f->constant, ctx->fqctx);
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
		success = _fq_nmod_mpoly_factor_separable(t, f->poly + j, ctx, sep);
		if (!success)
			goto cleanup;

        FLINT_ASSERT(fq_nmod_is_one(t->constant, ctx->fqctx));

        fq_nmod_mpoly_factor_fit_length(g, g->num + t->num, ctx);
        for (i = 0; i < t->num; i++)
        {
            fmpz_mul(g->exp + g->num, t->exp + i, f->exp + j);
            fq_nmod_mpoly_swap(g->poly + g->num, t->poly + i, ctx);
            g->num++;
        }
    }

    fq_nmod_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    fq_nmod_mpoly_factor_clear(t, ctx);
    fq_nmod_mpoly_factor_clear(g, ctx);

    FLINT_ASSERT(!success || fq_nmod_mpoly_factor_matches(A, f, ctx));

    return success;
}


int fq_nmod_mpoly_factor_squarefree(
    fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    return fq_nmod_mpoly_factor_separable(f, A, ctx, 0);
}

