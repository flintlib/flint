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
		t = nmod_pow_fmpz(t, e, ctx->ffinfo->mod);
		fac->constant = nmod_mul(fac->constant, t, ctx->ffinfo->mod);
		return;
	}
	else
	{
		nmod_mpoly_factor_append_fmpz(fac, a, e, ctx);
	}
}

static void nmod_mpoly_factor_mul_mpoly_ui(
	nmod_mpoly_factor_t fac,
	const nmod_mpoly_t a,
	ulong e,
	const nmod_mpoly_ctx_t ctx)
{
	if (nmod_mpoly_is_ui(a, ctx))
	{
		ulong t = nmod_mpoly_get_ui(a, ctx);
		t = nmod_pow_ui(t, e, ctx->ffinfo->mod);
		fac->constant = nmod_mul(fac->constant, t, ctx->ffinfo->mod);
		return;
	}
	else
	{
		nmod_mpoly_factor_append_ui(fac, a, e, ctx);
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

	a->constant = nmod_mul(b->constant, c->constant, ctx->ffinfo->mod);
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


static int _squarefree_factors(
	nmod_mpoly_factor_t f,
	const nmod_mpoly_t A,
	const nmod_mpoly_ctx_t ctx)
{
	int success;
	slong var;
    ulong k;
	fmpz * shift, * stride;
	nmod_mpoly_factor_t tempf;
	fmpz_t g, gr, p, pk;
    nmod_mpoly_t B, C, U, V, W, G;

    nmod_mpoly_init(B, ctx);
    nmod_mpoly_init(C, ctx);
    nmod_mpoly_init(U, ctx);
    nmod_mpoly_init(V, ctx);
    nmod_mpoly_init(W, ctx);
    nmod_mpoly_init(G, ctx);

	fmpz_init_set_ui(p, ctx->ffinfo->mod.n);
	fmpz_init(pk);
	fmpz_init(g);
	fmpz_init(gr);
	nmod_mpoly_factor_init(tempf, ctx);
	shift = _fmpz_vec_init(ctx->minfo->nvars);
	stride = _fmpz_vec_init(ctx->minfo->nvars);

	nmod_mpoly_factor_one(f, ctx);
	nmod_mpoly_set(C, A, ctx);

	for (var = 0; var < ctx->minfo->nvars; var++)
	{
        nmod_mpoly_derivative(G, C, var, ctx);
        success = nmod_mpoly_gcd_cofactors(C, W, V, C, G, ctx);
        if (!success)
            goto cleanup;

        for (k = 1; k + 1 < ctx->ffinfo->mod.n &&
                            !(nmod_mpoly_derivative(G, W, var, ctx),
                              nmod_mpoly_sub(U, V, G, ctx),
                              nmod_mpoly_is_zero(U, ctx)); k++)
        {
            success = nmod_mpoly_gcd_cofactors(G, W, V, W, U, ctx);
            if (!success)
	            goto cleanup;

            nmod_mpoly_factor_mul_mpoly_ui(f, G, k, ctx);

            if (!nmod_mpoly_is_one(W, ctx))
            {
                success = nmod_mpoly_divides(U, C, W, ctx);
                FLINT_ASSERT(success);
                nmod_mpoly_swap(C, U, ctx);
            }
        }

        nmod_mpoly_factor_mul_mpoly_ui(f, W, k, ctx);
	}

	if (nmod_mpoly_is_ui(C, ctx))
	{
		f->constant = nmod_mul(f->constant, nmod_mpoly_get_ui(C, ctx),
                                                             ctx->ffinfo->mod);
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
		success = _squarefree_factors(tempf, C, ctx);
		if (!success)
			goto cleanup;

        /* tempf = temp^pk */
	    tempf->constant = nmod_pow_fmpz(tempf->constant, pk, ctx->ffinfo->mod);
	    _fmpz_vec_scalar_mul_fmpz(tempf->exp, tempf->exp, tempf->num, pk);

        /* f = f * tempf */
		nmod_mpoly_factor_mul_pairwise_prime(f, f, tempf, ctx);
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
	nmod_mpoly_factor_clear(tempf, ctx);
	_fmpz_vec_clear(shift, ctx->minfo->nvars);
	_fmpz_vec_clear(stride, ctx->minfo->nvars);

	return success;
}



int nmod_mpoly_factor_squarefree(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, v;
    nmod_mpoly_t c;
    nmod_mpoly_univar_t u;
    nmod_mpoly_factor_t g, t;
    fmpz * var_powers;

    f->num = 0;
    if (nmod_mpoly_is_ui(A, ctx))
    {
		f->constant = nmod_mpoly_get_ui(A, ctx);
        return 1;
    }

    FLINT_ASSERT(A->length > 0);
	f->constant = A->coeffs[0];
	nmod_mpoly_factor_fit_length(f, 1, ctx);
	nmod_mpoly_make_monic(f->poly + 0, A, ctx);
	fmpz_one(f->exp + 0);
	f->num = 1;

    nmod_mpoly_factor_init(g, ctx);
    nmod_mpoly_factor_init(t, ctx);
    nmod_mpoly_univar_init(u, ctx);
    nmod_mpoly_init(c, ctx);
    var_powers = _fmpz_vec_init(ctx->minfo->nvars);

    FLINT_ASSERT(nmod_mpoly_factor_matches(A, f, ctx));

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        g->constant = f->constant;
        g->num = 0;
        for (j = 0; j < f->num; j++)
        {
            FLINT_ASSERT(fmpz_is_one(f->exp + j));

            nmod_mpoly_to_univar(u, f->poly + j, v, ctx);
            FLINT_ASSERT(u->length > 0);

            success = _nmod_mpoly_vec_content_mpoly(c, u->coeffs, u->length, ctx);
            if (!success)
                goto cleanup;

            fmpz_add(var_powers + v, var_powers + v, u->exps + u->length - 1);
            for (i = 0; i < u->length; i++)
            {
                fmpz_sub(u->exps + i, u->exps + i, u->exps + u->length - 1);
                success = nmod_mpoly_divides(u->coeffs + i, u->coeffs + i, c, ctx);
                FLINT_ASSERT(success);
            }

            nmod_mpoly_factor_mul_mpoly_ui(g, c, 1, ctx);

            if (u->length > 1)
            {
                nmod_mpoly_from_univar_bits(c, A->bits, u, v, ctx);
                nmod_mpoly_factor_append_ui(g, c, 1, ctx);
            }
            else
            {
                FLINT_ASSERT(nmod_mpoly_is_one(u->coeffs + 0, ctx));
            }
        }

        nmod_mpoly_factor_swap(f, g, ctx);
    }

	FLINT_ASSERT(nmod_mpoly_factor_is_pairwise_prime(f, ctx) != 0);

    g->constant = f->constant;
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
		success = _squarefree_factors(t, f->poly + j, ctx);
		if (!success)
			goto cleanup;

        FLINT_ASSERT(fmpz_is_one(f->exp + j));
        FLINT_ASSERT(1 == t->constant);

        nmod_mpoly_factor_fit_length(g, g->num + t->num, ctx);
        for (i = 0; i < t->num; i++)
        {
            fmpz_swap(g->exp + g->num, t->exp + i);
            nmod_mpoly_swap(g->poly + g->num, t->poly + i, ctx);
            g->num++;
        }
    }

    nmod_mpoly_factor_swap(f, g, ctx);

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        if (fmpz_is_zero(var_powers + v))
            continue;

        nmod_mpoly_factor_fit_length(f, f->num + 1, ctx);
        nmod_mpoly_gen(f->poly + f->num, v, ctx);
        fmpz_swap(f->exp + f->num, var_powers + v);
        f->num++;
    }

    success = 1;

cleanup:

    nmod_mpoly_factor_clear(t, ctx);
    nmod_mpoly_factor_clear(g, ctx);
    nmod_mpoly_univar_clear(u, ctx);
    nmod_mpoly_clear(c, ctx);
    _fmpz_vec_clear(var_powers, ctx->minfo->nvars);

    FLINT_ASSERT(!success || nmod_mpoly_factor_matches(A, f, ctx));

    return success;
}

