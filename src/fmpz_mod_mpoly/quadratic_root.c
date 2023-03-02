/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


int fmpz_mod_mpoly_quadratic_root(
    fmpz_mod_mpoly_t Q,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        fmpz_mod_mpoly_zero(Q, ctx);
        return 1;
    }

    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        return fmpz_mod_mpoly_sqrt(Q, B, ctx);
    }

    if (fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)))
    {
        nmod_mpoly_ctx_t nctx;
        nmod_mpoly_t nQ, nA, nB;

        nctx->minfo[0] = ctx->minfo[0];
        nmod_init(&nctx->mod, fmpz_get_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)));
        nmod_mpoly_init(nQ, nctx);
        nmod_mpoly_init(nA, nctx);
        nmod_mpoly_init(nB, nctx);

        _fmpz_mod_mpoly_get_nmod_mpoly(nA, nctx, A, ctx);
        _fmpz_mod_mpoly_get_nmod_mpoly(nB, nctx, B, ctx);
        success = nmod_mpoly_quadratic_root(nQ, nA, nB, nctx);
        _fmpz_mod_mpoly_set_nmod_mpoly(Q, ctx, nQ, nctx);

        nmod_mpoly_clear(nA, nctx);
        nmod_mpoly_clear(nQ, nctx);
        nmod_mpoly_clear(nB, nctx);

        return success;
    }
    else
    {
        fmpz_t c, c2;
        fmpz_mod_mpoly_t t1, t2;

        fmpz_init(c);
        fmpz_init(c2);

        FLINT_ASSERT(fmpz_is_odd(fmpz_mod_ctx_modulus(ctx->ffinfo)));

        fmpz_fdiv_q_2exp(c, fmpz_mod_ctx_modulus(ctx->ffinfo), 1);
        fmpz_mod_mul(c2, c, c, ctx->ffinfo);

        fmpz_mod_mpoly_init(t1, ctx);
        fmpz_mod_mpoly_init(t2, ctx);

        fmpz_mod_mpoly_mul(t1, A, A, ctx);
        fmpz_mod_mpoly_scalar_addmul_fmpz(t2, B, t1, c2, ctx);
        success = fmpz_mod_mpoly_sqrt(t1, t2, ctx);
        if (success)
            fmpz_mod_mpoly_scalar_addmul_fmpz(Q, t1, A, c, ctx);

        fmpz_mod_mpoly_clear(t1, ctx);
        fmpz_mod_mpoly_clear(t2, ctx);

        fmpz_clear(c);
        fmpz_clear(c2);

        return success;
    }
}
