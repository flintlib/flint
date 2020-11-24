/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


static void nmod_mpoly_factor_mul_mpoly_fmpz(
	nmod_mpoly_factor_t fac,
	const nmod_mpoly_t a,
	const fmpz_t e,
	const nmod_mpoly_ctx_t ctx)
{
	if (nmod_mpoly_is_ui(a, ctx))
	{
		ulong t = nmod_mpoly_get_ui(a, ctx);
		t = nmod_pow_fmpz(t, e, ctx->mod);
		fac->constant = nmod_mul(fac->constant, t, ctx->mod);
		return;
	}
	else
	{
		nmod_mpoly_factor_append_fmpz(fac, a, e, ctx);
	}
}


#if FLINT_WANT_ASSERT
/*
    return: 0 no
            1 yes
           -1 don't know
*/
static int nmod_mpoly_factor_is_pairwise_prime(
	const nmod_mpoly_factor_t f, 
	const nmod_mpoly_ctx_t ctx)
{
	int result;
	slong i, j;
	nmod_mpoly_t g;

	nmod_mpoly_init(g, ctx);

	for (i = 0; i + 1 < f->num; i++)
	for (j = i + 1; j < f->num; j++)
	{
		/* make sure factors are monic */
		if (f->poly[i].length < 1 ||
            f->poly[j].length < 1 ||
            f->poly[i].coeffs[0] != 1 ||
            f->poly[j].coeffs[0] != 1)
		{
			result = 0;
			goto cleanup;
		}

		if (nmod_mpoly_gcd(g, f->poly + i, f->poly + j, ctx))
		{
			if (!nmod_mpoly_is_one(g, ctx))
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

	nmod_mpoly_clear(g, ctx);

	return result;
}
#endif

/*
	b and c are pairwise prime
	produce a=b*c pairwise prime and return 1
		else return 0 with a undefined
*/
static int nmod_mpoly_factor_mul_pairwise_prime(
	nmod_mpoly_factor_t a,
	const nmod_mpoly_factor_t b,
	const nmod_mpoly_factor_t c,
	const nmod_mpoly_ctx_t ctx)
{
	int success;
	slong i, j;
	nmod_mpoly_t T1, T2;
	nmod_mpoly_struct * g;
	fmpz_t t;

	if (a == b || a == c)
	{
		nmod_mpoly_factor_t ta;
		nmod_mpoly_factor_init(ta, ctx);
		success = nmod_mpoly_factor_mul_pairwise_prime(ta, b, c, ctx);
		nmod_mpoly_factor_swap(a, ta, ctx);
		nmod_mpoly_factor_clear(ta, ctx);
		return success;
	}

	fmpz_init(t);
	nmod_mpoly_init(T1, ctx);
	nmod_mpoly_init(T2, ctx);

	FLINT_ASSERT(nmod_mpoly_factor_is_pairwise_prime(b, ctx) == 1);
	FLINT_ASSERT(nmod_mpoly_factor_is_pairwise_prime(c, ctx) == 1);

	g = (nmod_mpoly_struct *) flint_malloc(b->num*c->num*sizeof(nmod_mpoly_struct));
	/* g[i,j] = gcd(b[i], c[j]) */
	for (i = 0; i < b->num; i++)
	for (j = 0; j < c->num; j++)
		nmod_mpoly_init(g + i*c->num + j, ctx);

	a->constant = nmod_mul(b->constant, c->constant, ctx->mod);
	a->num = 0;

	for (i = 0; i < b->num; i++)
	for (j = 0; j < c->num; j++)
	{
		if (!nmod_mpoly_gcd(g + i*c->num + j, b->poly + i, c->poly + j, ctx))
		{
			success = 0;
			goto cleanup;
		}

		fmpz_add(t, b->exp + i, c->exp + j);
		nmod_mpoly_factor_mul_mpoly_fmpz(a, g + i*c->num + j, t, ctx);
	}

	for (i = 0; i < b->num; i++)
	{
		nmod_mpoly_set(T1, b->poly + i, ctx);
		for (j = 0; j < c->num; j++)
		{
			success = nmod_mpoly_divides(T1, T1, g + i*c->num + j, ctx);
            FLINT_ASSERT(success);
		}
		nmod_mpoly_factor_mul_mpoly_fmpz(a, T1, b->exp + i, ctx);
	}	

	for (j = 0; j < c->num; j++)
	{
		nmod_mpoly_set(T1, c->poly + j, ctx);
		for (i = 0; i < b->num; i++)
		{
			success = nmod_mpoly_divides(T1, T1, g + i*c->num + j, ctx);
            FLINT_ASSERT(success);
		}
		nmod_mpoly_factor_mul_mpoly_fmpz(a, T1, c->exp + j, ctx);
	}	

	success = 1;

cleanup:

	for (i = 0; i < b->num; i++)
	for (j = 0; j < c->num; j++)
		nmod_mpoly_clear(g + i*c->num + j, ctx);

	flint_free(g);

	nmod_mpoly_clear(T1, ctx);
	nmod_mpoly_clear(T2, ctx);
	fmpz_clear(t);

	if (success)
	{
		nmod_mpoly_t ae, be, ce;
		nmod_mpoly_init(ae, ctx);
		nmod_mpoly_init(be, ctx);
		nmod_mpoly_init(ce, ctx);
		FLINT_ASSERT(nmod_mpoly_factor_is_pairwise_prime(a, ctx) == 1);
		nmod_mpoly_factor_expand(be, b, ctx);
		nmod_mpoly_factor_expand(ce, c, ctx);
		nmod_mpoly_mul(ae, be, ce, ctx);
		FLINT_ASSERT(nmod_mpoly_factor_matches(ae, a, ctx));
		nmod_mpoly_clear(ae, ctx);
		nmod_mpoly_clear(be, ctx);
		nmod_mpoly_clear(ce, ctx);
	}
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
    nmod_mpoly_factor_t f,
    nmod_mpoly_t a,
    ulong k,
    int * vars_left,
    const nmod_mpoly_ctx_t ctx,
    int sep,
    nmod_mpoly_t t) /* temp */
{
    slong v, org = f->num;

    if (nmod_mpoly_is_ui(a, ctx))
    {
        FLINT_ASSERT(nmod_mpoly_is_one(a, ctx));
        return 1;
    }

    nmod_mpoly_factor_fit_length(f, org + 1, ctx);
    nmod_mpoly_swap(f->poly + org, a, ctx);
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
            nmod_mpoly_derivative(t, f->poly + i, v, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
            {
                /* f[i] has zero derivative */
                i++;
                continue;
            }

            nmod_mpoly_factor_fit_length(f, f->num + 1, ctx);
            fmpz_set_ui(f->exp + f->num, k);

            if (!nmod_mpoly_gcd_cofactors(f->poly + f->num, f->poly + i, t,
                                                          f->poly + i, t, ctx))
            {
                return 0;
            }

            if (nmod_mpoly_is_ui(f->poly + f->num, ctx))
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


int _nmod_mpoly_factor_separable(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx,
    int sep)
{
	int success;
	slong v, var, j;
    ulong k;
    int * vars_left;
	fmpz * shift, * stride;
	nmod_mpoly_factor_t Tf;
	fmpz_t g, gr, p, pk;
    nmod_mpoly_t B, C, U, V, W, G;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0] == 1);

    nmod_mpoly_init(B, ctx);
    nmod_mpoly_init(C, ctx);
    nmod_mpoly_init(U, ctx);
    nmod_mpoly_init(V, ctx);
    nmod_mpoly_init(W, ctx);
    nmod_mpoly_init(G, ctx);

	fmpz_init_set_ui(p, ctx->mod.n);
	fmpz_init(pk);
	fmpz_init(g);
	fmpz_init(gr);
	nmod_mpoly_factor_init(Tf, ctx);
	shift = _fmpz_vec_init(ctx->minfo->nvars);
	stride = _fmpz_vec_init(ctx->minfo->nvars);
    vars_left = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, int);
    for (v = 0; v < ctx->minfo->nvars; v++)
        vars_left[v] = 1;

    nmod_mpoly_factor_one(f, ctx);
    nmod_mpoly_set(C, A, ctx);

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        /* take the next variable with shortest derivative */
        var = -1;
        for (v = 0; v < ctx->minfo->nvars; v++)
        {
            if (!vars_left[v])
                continue;

            nmod_mpoly_derivative(U, C, v, ctx);

            if (var < 0 || U->length < G->length)
            {
                var = v;
                nmod_mpoly_swap(G, U, ctx);
            }
        }

        FLINT_ASSERT(var >= 0);
        FLINT_ASSERT(vars_left[var] == 1);

        vars_left[var] = 0;

        success = nmod_mpoly_gcd_cofactors(C, W, V, C, G, ctx);
        if (!success)
            goto cleanup;

        for (k = 1; k + 1 < ctx->mod.n &&
                            !(nmod_mpoly_derivative(G, W, var, ctx),
                              nmod_mpoly_sub(U, V, G, ctx),
                              nmod_mpoly_is_zero(U, ctx)); k++)
        {
            success = nmod_mpoly_gcd_cofactors(G, W, V, W, U, ctx);
            if (!success)
	            goto cleanup;

            success = _append_factor_sep(f, G, k, vars_left, ctx, sep, U);
            if (!success)
                goto cleanup;

            if (!nmod_mpoly_is_one(W, ctx))
            {
                success = nmod_mpoly_divides(U, C, W, ctx);
                FLINT_ASSERT(success);
                nmod_mpoly_swap(C, U, ctx);
            }
        }

        success = _append_factor_sep(f, W, k, vars_left, ctx, sep, U);
        if (!success)
            goto cleanup;
	}

    if (nmod_mpoly_is_ui(C, ctx))
	{
        FLINT_ASSERT(nmod_mpoly_is_one(C, ctx));
	}
	else
	{
		nmod_mpoly_deflation(shift, stride, C, ctx);
		fmpz_zero(g);
		for (var = 0; var < ctx->minfo->nvars; var++)
		{
			fmpz_gcd(g, g, stride + var);
			fmpz_gcd(g, g, shift + var);
		}

		fmpz_pow_ui(pk, p, fmpz_remove(gr, g, p));
        FLINT_ASSERT(fmpz_cmp_ui(pk, 1) > 0);

		for (var = 0; var < ctx->minfo->nvars; var++)
		{
			fmpz_set(stride + var, pk);
			fmpz_zero(shift + var);
		}

		nmod_mpoly_deflate(C, C, shift, stride, ctx);
		success = _nmod_mpoly_factor_separable(Tf, C, ctx, sep);
		if (!success)
			goto cleanup;

        /* f *= Tf^pk */
        FLINT_ASSERT(Tf->constant == 1);
	    _fmpz_vec_scalar_mul_fmpz(Tf->exp, Tf->exp, Tf->num, pk);
		nmod_mpoly_factor_mul_pairwise_prime(f, f, Tf, ctx);
	}

	success = 1;

cleanup:

    nmod_mpoly_clear(C, ctx);
    nmod_mpoly_clear(U, ctx);
    nmod_mpoly_clear(V, ctx);
    nmod_mpoly_clear(W, ctx);
    nmod_mpoly_clear(G, ctx);

	fmpz_clear(p);
	fmpz_clear(pk);
	fmpz_clear(g);
	fmpz_clear(gr);
	nmod_mpoly_factor_clear(Tf, ctx);
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
int nmod_mpoly_factor_separable(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx,
    int sep)
{
    int success;
    slong i, j;
    nmod_mpoly_factor_t g, t;

    if (!nmod_mpoly_factor_content(f, A, ctx))
        return 0;

    nmod_mpoly_factor_init(g, ctx);
    nmod_mpoly_factor_init(t, ctx);

    g->constant = f->constant;
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
		success = _nmod_mpoly_factor_separable(t, f->poly + j, ctx, sep);
		if (!success)
			goto cleanup;

        FLINT_ASSERT(1 == t->constant);

        nmod_mpoly_factor_fit_length(g, g->num + t->num, ctx);
        for (i = 0; i < t->num; i++)
        {
            fmpz_mul(g->exp + g->num, t->exp + i, f->exp + j);
            nmod_mpoly_swap(g->poly + g->num, t->poly + i, ctx);
            g->num++;
        }
    }

    nmod_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    nmod_mpoly_factor_clear(t, ctx);
    nmod_mpoly_factor_clear(g, ctx);

    FLINT_ASSERT(!success || nmod_mpoly_factor_matches(A, f, ctx));

    return success;
}


int nmod_mpoly_factor_squarefree(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_factor_separable(f, A, ctx, 0);
}

