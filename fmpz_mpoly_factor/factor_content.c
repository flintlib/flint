/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"

/*
    return factors that are primitive wrt each variable
*/
int fmpz_mpoly_factor_content(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong j, v;
    fmpz_mpoly_t c;
    fmpz_mpoly_univar_t u;
    fmpz_mpoly_factor_t g;
    fmpz * var_powers;

    f->num = 0;
    if (fmpz_mpoly_is_fmpz(A, ctx))
    {
        fmpz_mpoly_get_fmpz(f->constant, A, ctx);
        return 1;
    }

    FLINT_ASSERT(A->length > 0);
    _fmpz_vec_content(f->constant, A->coeffs, A->length);
    if (fmpz_sgn(A->coeffs + 0) < 0)
        fmpz_neg(f->constant, f->constant);

    fmpz_mpoly_factor_fit_length(f, 1, ctx);
    fmpz_mpoly_scalar_divexact_fmpz(f->poly + 0, A, f->constant, ctx);
    fmpz_one(f->exp + 0);
    f->num = 1;

    fmpz_mpoly_factor_init(g, ctx);
    fmpz_mpoly_univar_init(u, ctx);
    fmpz_mpoly_init(c, ctx);
    var_powers = _fmpz_vec_init(ctx->minfo->nvars);

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        fmpz_swap(g->constant, f->constant);
        g->num = 0;
        for (j = 0; j < f->num; j++)
        {
            FLINT_ASSERT(fmpz_is_one(f->exp + j));

            fmpz_mpoly_to_univar(u, f->poly + j, v, ctx);

            fmpz_add(var_powers + v, var_powers + v, u->exps + u->length - 1);
            _mpoly_gen_shift_right_fmpz(f->poly[j].exps, f->poly[j].bits,
                    f->poly[j].length, v, u->exps + u->length - 1, ctx->minfo);

            if (u->length < 2)
            {
                FLINT_ASSERT(u->length == 1);
                fmpz_mpoly_factor_fit_length(g, g->num + 1, ctx);
                fmpz_swap(g->exp + g->num, f->exp + j);
                fmpz_mpoly_swap(g->poly + g->num, f->poly + j, ctx);
                g->num += !fmpz_mpoly_is_fmpz(g->poly + g->num, ctx);
                continue;
            }

            success = _fmpz_mpoly_vec_content_mpoly(c, u->coeffs, u->length, ctx);
            if (!success)
                goto cleanup;

            if (fmpz_mpoly_is_fmpz(c, ctx))
            {
                fmpz_mpoly_factor_fit_length(g, g->num + 1, ctx);
                fmpz_swap(g->exp + g->num, f->exp + j);
                fmpz_mpoly_swap(g->poly + g->num, f->poly + j, ctx);
                g->num++;
            }
            else
            {
                fmpz_mpoly_factor_fit_length(g, g->num + 2, ctx);
                success = fmpz_mpoly_divides(g->poly + g->num, f->poly + j, c, ctx);
                FLINT_ASSERT(success);
                fmpz_mpoly_swap(g->poly + g->num + 1, c, ctx);
                fmpz_one(g->exp + g->num);
                fmpz_one(g->exp + g->num + 1);
                g->num += 2;
            }
        }

        fmpz_mpoly_factor_swap(f, g, ctx);
    }

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        if (fmpz_is_zero(var_powers + v))
            continue;

        fmpz_mpoly_factor_fit_length(f, f->num + 1, ctx);
        fmpz_mpoly_gen(f->poly + f->num, v, ctx);
        fmpz_swap(f->exp + f->num, var_powers + v);
        f->num++;
    }

    success = 1;

cleanup:

    fmpz_mpoly_factor_clear(g, ctx);
    fmpz_mpoly_univar_clear(u, ctx);
    fmpz_mpoly_clear(c, ctx);
    _fmpz_vec_clear(var_powers, ctx->minfo->nvars);

    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));

    return success;
}

