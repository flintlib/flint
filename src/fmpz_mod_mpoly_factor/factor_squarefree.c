/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"
#include "nmod_mpoly_factor.h"

static void _fmpz_mod_mpoly_factor_mul_mpoly_fmpz(
    fmpz_mod_mpoly_factor_t f,
    fmpz_mod_mpoly_t A,
    const fmpz_t e,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_fmpz(A, ctx))
    {
        FLINT_ASSERT(fmpz_mod_mpoly_is_one(A, ctx));
    }
    else
    {
        fmpz_mod_mpoly_factor_append_fmpz_swap(f, A, e, ctx);
    }
}

int _fmpz_mod_mpoly_factor_separable(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx,
    int sep)
{
	int success;
	slong v, var;
	fmpz_t k;
    fmpz_mod_mpoly_t U, V, W, G;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(fmpz_is_one(A->coeffs + 0));

    if (fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)))
    {
        nmod_mpoly_ctx_t nctx;
        nmod_mpoly_t nA;
        nmod_mpoly_factor_t nf;

        *nctx->minfo = *ctx->minfo;
        nmod_init(&nctx->mod, fmpz_get_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)));
        nmod_mpoly_init(nA, nctx);
        nmod_mpoly_factor_init(nf, nctx);

        _fmpz_mod_mpoly_get_nmod_mpoly(nA, nctx, A, ctx);
        success = _nmod_mpoly_factor_separable(nf, nA, nctx, sep);
        _fmpz_mod_mpoly_factor_set_nmod_mpoly_factor(f, ctx, nf, nctx);

        nmod_mpoly_factor_clear(nf, nctx);
        nmod_mpoly_clear(nA, nctx);

        return success;
    }

    fmpz_mod_mpoly_factor_one(f, ctx);

    if (!fmpz_mod_mpoly_degrees_fit_si(A, ctx))
        return 0;

    fmpz_init(k);
    fmpz_mod_mpoly_init(U, ctx);
    fmpz_mod_mpoly_init(V, ctx);
    fmpz_mod_mpoly_init(W, ctx);
    fmpz_mod_mpoly_init(G, ctx);

    /* take the variable with shortest derivative */
    var = -1;
    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        fmpz_mod_mpoly_derivative(U, A, v, ctx);

        if (U->length > 0 && (var < 0 || U->length < G->length))
        {
            var = v;
            fmpz_mod_mpoly_swap(V, U, ctx);
        }
    }

    if (var < 0)
    {
        FLINT_ASSERT(fmpz_mod_mpoly_is_one(A, ctx));
        success = 1;
        goto cleanup;
    }

    success = fmpz_mod_mpoly_gcd_cofactors(G, W, V, A, V, ctx);
    if (!success)
        goto cleanup;

    fmpz_one(k);
 
    while (1)
    {
        fmpz_add_ui(k, k, 1);
        if (fmpz_cmp(k, fmpz_mod_ctx_modulus(ctx->ffinfo)) >= 0)
            break;
        fmpz_sub_ui(k, k, 1);

        fmpz_mod_mpoly_derivative(G, W, var, ctx);
        fmpz_mod_mpoly_sub(U, V, G, ctx);
        if (fmpz_mod_mpoly_is_zero(U, ctx))
            break;

        success = fmpz_mod_mpoly_gcd_cofactors(G, W, V, W, U, ctx);
        if (!success)
            goto cleanup;

        _fmpz_mod_mpoly_factor_mul_mpoly_fmpz(f, G, k, ctx);

        fmpz_add_ui(k, k, 1);
    }

    _fmpz_mod_mpoly_factor_mul_mpoly_fmpz(f, W, k, ctx);

	success = 1;

cleanup:

    fmpz_clear(k);
    fmpz_mod_mpoly_clear(U, ctx);
    fmpz_mod_mpoly_clear(V, ctx);
    fmpz_mod_mpoly_clear(W, ctx);
    fmpz_mod_mpoly_clear(G, ctx);

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
int fmpz_mod_mpoly_factor_separable(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx,
    int sep)
{
    int success;
    slong i, j;
    fmpz_mod_mpoly_factor_t g, t;

    if (!fmpz_mod_mpoly_factor_content(f, A, ctx))
        return 0;

    fmpz_mod_mpoly_factor_init(g, ctx);
    fmpz_mod_mpoly_factor_init(t, ctx);

    fmpz_swap(g->constant, f->constant);
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
		success = _fmpz_mod_mpoly_factor_separable(t, f->poly + j, ctx, sep);
		if (!success)
			goto cleanup;

        FLINT_ASSERT(fmpz_is_one(t->constant));

        fmpz_mod_mpoly_factor_fit_length(g, g->num + t->num, ctx);
        for (i = 0; i < t->num; i++)
        {
            fmpz_mul(g->exp + g->num, t->exp + i, f->exp + j);
            fmpz_mod_mpoly_swap(g->poly + g->num, t->poly + i, ctx);
            g->num++;
        }
    }

    fmpz_mod_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    fmpz_mod_mpoly_factor_clear(t, ctx);
    fmpz_mod_mpoly_factor_clear(g, ctx);

    FLINT_ASSERT(!success || fmpz_mod_mpoly_factor_matches(A, f, ctx));

    return success;
}


int fmpz_mod_mpoly_factor_squarefree(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_factor_separable(f, A, ctx, 0);
}

