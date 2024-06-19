/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"
#include "gr.h"

void
fmpz_mpoly_q_print_pretty(const fmpz_mpoly_q_t f, const char ** x, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(f), ctx))
    {
        fmpz_mpoly_print_pretty(fmpz_mpoly_q_numref(f), x, ctx);
    }
    else if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(f), ctx))
    {
        flint_printf("(");
        fmpz_mpoly_print_pretty(fmpz_mpoly_q_numref(f), x, ctx);
        flint_printf(")/");
        fmpz_mpoly_print_pretty(fmpz_mpoly_q_denref(f), x, ctx);
    }
    else
    {
        flint_printf("(");
        fmpz_mpoly_print_pretty(fmpz_mpoly_q_numref(f), x, ctx);
        flint_printf(")/(");
        fmpz_mpoly_print_pretty(fmpz_mpoly_q_denref(f), x, ctx);
        flint_printf(")");
    }
}

char * fmpz_mpoly_q_get_str_pretty(const fmpz_mpoly_q_t f, const char ** vars, const fmpz_mpoly_ctx_t ctx)
{
    gr_ctx_t grctx;
    char * s;

    gr_ctx_init_fmpz_mpoly_q(grctx, ctx->minfo->nvars, ctx->minfo->ord);
    if (vars != NULL)
        GR_MUST_SUCCEED(gr_ctx_set_gen_names(grctx, vars));
    GR_MUST_SUCCEED(gr_get_str(&s, f, grctx));
    gr_ctx_clear(grctx);

    return s;
}

int
fmpz_mpoly_q_set_str_pretty(fmpz_mpoly_q_t res, const char * s, const char ** vars, fmpz_mpoly_ctx_t ctx)
{
    gr_ctx_t grctx;
    int ret;

    gr_ctx_init_fmpz_mpoly_q(grctx, ctx->minfo->nvars, ctx->minfo->ord);
    if (vars != NULL)
        GR_MUST_SUCCEED(gr_ctx_set_gen_names(grctx, vars));
    ret = (GR_SUCCESS == gr_set_str(res, s, grctx)) ? 0 : -1;

    gr_ctx_clear(grctx);
    return ret;
}
