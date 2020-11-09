/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

/* return factors that are primitive wrt each variable */
int nmod_mpoly_factor_content(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong nvars = ctx->minfo->nvars;
    slong i, j, k, v;
    nmod_mpoly_t c;
    nmod_mpoly_univar_t u;
    nmod_mpoly_factor_t g;
    fmpz * pows;
    slong * perm;

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
    nmod_mpoly_univar_init(u, ctx);
    nmod_mpoly_init(c, ctx);
    pows = _fmpz_vec_init(nvars);

    perm = FLINT_ARRAY_ALLOC(nvars, slong);
    for (k = nvars - 1; k >= 0; k--)
        perm[k] = k;

    mpoly_degrees_ffmpz(pows, A->exps, A->length, A->bits, ctx->minfo);

    for (i = 1; i < nvars; i++)
        for (j = i; j > 0 && fmpz_cmp(pows + perm[j], pows + perm[j-1]) < 0; j--)
            SLONG_SWAP(perm[j], perm[j - 1]);

    _fmpz_vec_zero(pows, nvars);

    for (k = ctx->minfo->nvars - 1; k >= 0; k--)
    {
        v = perm[k];

        g->constant = f->constant;
        g->num = 0;
        for (j = 0; j < f->num; j++)
        {
            FLINT_ASSERT(fmpz_is_one(f->exp + j));

            nmod_mpoly_to_univar(u, f->poly + j, v, ctx);
            FLINT_ASSERT(u->length > 0);

            fmpz_add(pows + v, pows + v, u->exps + u->length - 1);
            _mpoly_gen_shift_right_fmpz(f->poly[j].exps, f->poly[j].bits,
                    f->poly[j].length, v, u->exps + u->length - 1, ctx->minfo);

            if (u->length < 2)
            {
                FLINT_ASSERT(u->length == 1);
                nmod_mpoly_factor_fit_length(g, g->num + 1, ctx);
                fmpz_swap(g->exp + g->num, f->exp + j);
                nmod_mpoly_swap(g->poly + g->num, f->poly + j, ctx);
                g->num += !nmod_mpoly_is_ui(g->poly + g->num, ctx);
                continue;
            }

            success = _nmod_mpoly_vec_content_mpoly(c, u->coeffs, u->length, ctx);
            if (!success)
                goto cleanup;

            if (nmod_mpoly_is_ui(c, ctx))
            {
                nmod_mpoly_factor_fit_length(g, g->num + 1, ctx);
                fmpz_swap(g->exp + g->num, f->exp + j);
                nmod_mpoly_swap(g->poly + g->num, f->poly + j, ctx);
                g->num++;
            }
            else
            {
                nmod_mpoly_factor_fit_length(g, g->num + 2, ctx);
                success = nmod_mpoly_divides(g->poly + g->num, f->poly + j, c, ctx);
                FLINT_ASSERT(success);
                nmod_mpoly_swap(g->poly + g->num + 1, c, ctx);
                fmpz_one(g->exp + g->num);
                fmpz_one(g->exp + g->num + 1);
                g->num += 2;
            }
        }

        nmod_mpoly_factor_swap(f, g, ctx);
    }

    for (v = 0; v < nvars; v++)
    {
        if (fmpz_is_zero(pows + v))
            continue;

        nmod_mpoly_factor_fit_length(f, f->num + 1, ctx);
        nmod_mpoly_gen(f->poly + f->num, v, ctx);
        fmpz_swap(f->exp + f->num, pows + v);
        f->num++;
    }

    success = 1;

cleanup:

    nmod_mpoly_factor_clear(g, ctx);
    nmod_mpoly_univar_clear(u, ctx);
    nmod_mpoly_clear(c, ctx);
    _fmpz_vec_clear(pows, nvars);
    flint_free(perm);

    FLINT_ASSERT(!success || nmod_mpoly_factor_matches(A, f, ctx));

    return success;
}

